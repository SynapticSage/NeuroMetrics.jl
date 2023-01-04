
using DimensionalData
using ProgressMeter
using Plots
include(scriptsdir("causality", "init_trig_plot.jl"))

# Pull out the grid
grd = storage["grd"];
props = grd.props
dims = Symbol.((props..., "N", "horizon"))

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
        thingfunc=identity
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
            #@infiltrate
        catch
            @infiltrate
        end
    end

    P
end

thingfunc = x->mean((!).(isnan.(x)),dims=(:startWell,:stopWell,:hatrajnum))
P = block_of_plots(ca1pfc; 
                   thingfunc, 
                   plot_dims=(:cuemem, :correct))


thingfunc = x->sum(x, 
                    dims=(:startWell,:stopWell,:hatrajnum))
P = block_of_plots(cnts; thingfunc, 
                   plot_dims=(:cuemem, :correct))
plot(P..., layout=grid(3,3), titlefontsize=4)


Ps = [( 
  cntsub = cnts[startWell=i, stopWell=j];
    thingfunc = x->sum(x, 
                        dims=(:hatrajnum, :N));
    P = block_of_plots(cntsub; thingfunc, 
                       plot_dims=());
    plot(P..., layout=grid(3,3), titlefontsize=4)
   ) for (i,j) in Iterators.product(1:6, 1:6)]
plot(Ps...)
