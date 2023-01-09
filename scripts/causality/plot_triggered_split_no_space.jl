
using DimensionalData
using ProgressMeter
using Plots
using LazyGrids
import Plot
import Utils: namedtup
include(scriptsdir("causality", "init_trig_plot.jl"))

Plot.setparentfolder(plotsdir("causality", "plot_triggered_split_no_space.jl"))

# Pull out the grid
grd = storage["grd"];
props = grd.props
dims = Symbol.((props..., "N", "horizon"))
Plot.setappend(namedtup.tostring((;props)))

# Pull out checkpoint
checkpoint, cnts = storage["checkpoint"], storage["counts"];
countmap(values(checkpoint))

ca1pfc, pfcca1, cnts = DimArray(ca1pfc, dims), 
                 DimArray(pfcca1, dims),
                 DimArray(cnts, dims[1:length(props)+1]);
ca1pfc

"""
    block_of_plots

Return a block of plots over dimensions `plot_dims` from DimArray `thing`
"""
function block_of_plots(thing; 
        plot_dims=(:cuemem, :correct, :hatrajnum,:startWell,:stopWell),
        thingfunc=identity,
        labelnums=true, textlabelsize=4,
    )

    thing = thingfunc(thing)
    dims = name(thing.dims)
    sizes = Dict(dim=>size(thing,dim) for dim in dims)
    ranges = OrderedDict(dim=>1:size(thing,dim) for dim in dims)

    init = Array{Union{Plots.Plot,Missing}, 
                 length(plot_dims)}(missing,[sizes[d] for d in plot_dims]...);
    P = DimArray(init, plot_dims)

    @showprogress for vals in Iterators.product(collect.(getindex.([ranges], plot_dims))...)
        index = (;zip(plot_dims,vals)...)
        try
            T = !isempty(index) ? Utils.squeeze(getindex(thing; index...)) :
                                Utils.squeeze(thing)
            P[vals...] = plot(T; seriestype=(
                                  if ndims(T) == 1
                                      :line
                                  elseif ndims(T) == 2
                                      :heatmap
                                  end)
                             )
            if labelnums
                if ndims(T) == 1
                    coordsx = 1:length(thing)
                    A = [(x,y,t) for (x,y,t) in 
                         zip(coordsx[:], thing[:], text.(string.(thing[:]), textlabelsize))]
                    annotate!(A)
                elseif ndims(T) == 2
                    coordsy, coordsx = ndgrid(UnitRange.([1], size(thing))...)
                    A = [(x,y,t) for (x,y,t) in 
                         zip(coordsx[:], coordsy[:], text.(string.(thing[:]),textlabelsize))]
                    annotate!(A)
                end
            end
            #@infiltrate
        catch
            @infiltrate
        end
    end

    P
end

thingfunc = x->mean((!).(isnan.(x)), dims=(:startWell,:stopWell,:hatrajnum))
P = block_of_plots(ca1pfc; 
                   thingfunc, 
                   plot_dims=(:cuemem, :correct))

thingfunc = x->sum(x, 
                    dims=(:startWell,:stopWell,:hatrajnum))
P = block_of_plots(cnts; thingfunc, 
                   plot_dims=(:cuemem, :correct))
plot(P..., layout=grid(3,3), titlefontsize=4)


# Visualize counts per start and stop well
Ps = [( 
  cntsub = cnts[startWell=i, stopWell=j];
    thingfunc = x->sum(x, 
                        dims=(:hatrajnum, :N));
    P = block_of_plots(cntsub; thingfunc, 
                       plot_dims=());
    plot(P..., layout=grid(3,3), titlefontsize=4)
   ) for (i,j) in Iterators.product(1:6, 1:6)]
plot(Ps...; xticks=[], yticks=[], titlefontsize=2, colorbar=nothing, labelfontsize=2)
Plot.save("counts_by_startWell-stopWell")


