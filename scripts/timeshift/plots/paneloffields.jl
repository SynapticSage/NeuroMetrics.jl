    using Timeshift
    using Timeshift.types
    using Timeshift.shiftmetrics
    using Field.metrics
    using Plot
    using Plot.receptivefield
    using Utils.namedtup
    using DimensionalData
    using ProgressMeter
    import Plot
    using JLD2
    getshiftda(M::DimArray, s::Real) = M[:, M.dims[2].==s];
    function getshiftda(M::DimArray, s::Symbol)
        B = [findfirst(row) for row in eachrow(unitshift[:,1][:bestshift_bitsperspike].data .== shifts')]
        shiftB = [unit[b] for (unit,b) in zip(eachrow(unitshift),B)]
    end
    
    animal, day = "RY16", 36
    #animal, day = "RY22", 21
    #D=JLD2.load("/home/ryoung/tmp.jld2")
    #Utils.dict.load_dict_to_module!(Main,D)


    Plot.setfolder("cells","panel_of_fields")
    Plot.setappend("$animal-$day")

@time begin
    Plot.setappend("$animal-$day")

    @time F = load_fields();
    @time f = F[bestpartialmatch(keys(F), (;datacut=:all, widths=5, animal, day))];
    @time M = Timeshift.types.matrixform(f);
    @time m = M[:, M.dims[2].==0];
    size(m)
    mm = m[1:10]
    @time spikes, beh, ripples, cells = Load.load(animal,day);
    shifts = getshifts(f)
    nothing
end

# Setup basic metrics, like the area and dimensions of the dimarray
push_metric!.([M], [Field.metrics.meanrate, Field.metrics.maxrate, Field.metrics.bitsperspike, Field.metrics.coherence, Field.metrics.centroid])
push_celltable!( M, cells, :unit, :area)
push_shiftmetric!(M, best_tau!; metric=:bitsperspike)
push_dims!(M)
push_dims!(M, vec(getindex.(m, :coherence)); dim=:unit, metric=:coh_at_zero)

# Prefilter our fields : we could do this in the plot, but, we want a stable set over these figures
metricfilter = x->
(x[:area] .=="CA1" .&& x[:meanrate] .> 0.01 && x[:meanrate] .< 7 .&& x[:maxrate] < 35) .||
(x[:area] .=="PFC" .&& x[:meanrate] .> 0.005 .&& x[:maxrate] < 100)
area_conditional_c = Dict("CA1"=>:linear_kryw_5_100_c64_n256,
                          "PFC"=>:linear_kbc_5_95_c73_n256)
inds = vec(any(hcat(
                    [[metricfilter(mm) for mm in m]
                     for m in eachslice(M, dims=:unit)]...), dims=1))
Msub = M[inds,:]
Msub0 = getshiftda(Msub, 0)
@time push_metric!(Msub0, metrics.coherence;               name=:coherence, prog=true, skip_edge=true)

# Push different versions of my metric
@time push_metric!(Msub0, metrics.bcoherence;              name=:bcoherence, prog=true, skip_edge=true)
@time push_metric!(Msub0, metrics.blake_coherence;         name=:blake_coherence, prog=true)
@time push_metric!(Msub0, metrics.jake_coherence_pearson;  name=:jake_coherence_pearson, prog=true)
@time push_metric!(Msub0, metrics.jake_coherence_spearman; name=:jake_coherence_spearman, prog=true)


metriccompare = false
if metriccompare
    plts=[]; 
    fs=[:jake_coherence_spearman, :jake_coherence_pearson, :blake_coherence, :coherence, :bcoherence];
    theta = 0
    for (x,y) in Iterators.product(fs,fs); 
    push!(plts,
          scatter(Msub0[x], Msub0[y], camera=(theta,30), xlabel=string(x), ylabel=string(y))
         ); 
    plot!(0:2,0:2, c=:white,linestyle=:dash,linewidth=3, ylim=(0,2), xlim=(0,2))
    end
    plot(plts..., size=(1200,1200))
end

#push_metric!(Msub0, metrics.jake_coherence;  name=:jake_coherence2, prog=true)
push_celltable!(Msub0, cells, :unit, :area)
push_shiftmetric!(Msub0, best_tau!; metric=:bitsperspike)

#;pushover-cli "Readying to plot"
#savefile="/home/ryoung/tmp-$animal-$day.jld2"
#JLD2.@save(savefile)

# Display the fields at a zero shift
# TODO CAUTION can and will crash from field size!
    Pfieldsplay = fieldsplay(Msub0[1:60]; colorbar=false, 
                    metricdisp=[:unit, :area, :meanrate, :maxrate],
                    scale_spg=2, metricfilter, totalarea_aspect=1,
                    totalcount=100,
                    srt=[:area, :coherence], transpose=false,
                    title_width=20,area_conditional_c)
#serialize(expanduser("~/tmp.jld2"), (;Pfieldsplay, Msub0, metricfilter, area_conditional_c))
Plot.save("fields_at_zero_shift_srt=coherence")

Plot.setfolder("timeshifts", "field x shift", "$animal-$day")
for (u, shift) in Iterators.product(1:size(M,1),1:size(M,2))
    area = M[u, shift][:area]
    unit = M[u, shift][:unit]
    if area == "CA1"
        continue
    end
    @info "cell" area

    c =  area_conditional_c[area]
    plot(M[u,shift]; colorbar=false, 
                metricdisp=[:unit, :area, :meanrate, :maxrate],
                transpose=false, c)
    Plot.save((;unit,shift))
end

# GET A GIF OF ACTIVITY
prog = Progress(length(shifts);
                desc="Shifts")
anim = @animate for shift in shifts

    m = vec(getshiftda(Msub, shift))

    fieldsplay(m; colorbar=false, 
               metricdisp=[:area, :coherence, :meanrate, :coh_at_zero], 
               scale_spg=2, 
               totalarea_aspect=1,
               srt=[:area,:coh_at_zero],
               title_width=20,
               area_conditional_c)

    next!(prog)
end
gif(anim)

push_metric!(Msub0, Field.metrics.convexhull)

function plothull(U)
    U = U[:convexhull]
    hull = [VPolygon(U[:hullseg_grid][i])
        for i in 1:length(U[:hullseg_grid])][Utils.na, :]
    if lenth(hull)>=1
    plot!(hull[1])
    else
    plot!()
    end
end
for m in Msub0
    plothull(Msub0)
end

# GET FRAME AT BEST
push_shiftmetric!(Msub, best_tauind!; metric=:bitsperspike)
bestshifts = (M[:,1])[:bestshiftind_bitsperspike]
m = [M[row, bestshifts[row]] for row in 1:size(M,1)]

fieldsplay(m; colorbar=false, 
           metricdisp=[:area, :shift,
                       :meanrate, :coh_at_zero], 
           scale_spg=2, 
           totalarea_aspect=1,
           metricfilter,
           srt=[:area,:coh_at_zero],
           title_width=20,
           area_conditional_c)
Plot.save("fields_at_BEST_shift_srt=coh_at_zero")

Plot.setfolder("timeshift","population","panel_of_fields", "best_sorted_differently")
fieldsplay(m; colorbar=false, 
           metricdisp=[:area, :shift,
                       :meanrate, :coh_at_zero], 
           scale_spg=2, 
           totalarea_aspect=1,
           metricfilter,
           srt=[:area,:coherence],
           title_width=20,
           area_conditional_c)
Plot.save("fields_at_BEST_shift_srt=>coherence")


fieldsplay(m; colorbar=false, 
           metricdisp=[:area,     :shift,
                       :meanrate, :coh_at_zero], 
           scale_spg=2, 
           totalarea_aspect=1,
           metricfilter,
           srt=[:area, :bitsperspike],
           title_width=20,
           area_conditional_c)

Plot.save("fields_at_BEST_shift_srt=>bitsperspike")
