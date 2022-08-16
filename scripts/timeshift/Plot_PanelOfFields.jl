
@time begin
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
    Plot.setfolder("timeshift","population","panel_of_fields")

    @time F = load_fields();
    @time f = F[bestpartialmatch(keys(F), (;datacut=:all, widths=5))];
    @time f = ShiftedFields(f);
    @time M = Timeshift.types.matrixform(f);
    @time m = M[:, M.dims[2].==0];
    size(m)
    mm = m[1:10]
    @time spikes, beh, ripples, cells = Load.load("RY16", 36);
    getshift(M::DimArray, s) = M[:, M.dims[2].==s];
    shifts = getshifts(f)
    nothing
end

# Setup basic metrics, like the area and dimensions of the dimarray
push_celltable!( M, cells, :unit, :area)
push_shiftmetric!(M, best_tau!; metric=:bitsperspike)
push_dims!(M); dim=:unit, metric=:coh_at_zero)
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

# Display the fields at a zero shift
Plot.setfolder("timeshift","population","panel_of_fields")
fieldsplay(getshift(Msub, 0); colorbar=false, metricdisp=[:area, :coherence, :meanrate,:unit],
           scale_spg=2, metricfilter, totalarea_aspect=1,
           srt=[:area,:coh_at_zero],
           title_width=20,area_conditional_c)
Plot.save("fields_at_zero_shift_srt=coherence_at_zero")


# GET A GIF OF ACTIVITY
prog = Progress(length(shifts);
                desc="Shifts")
anim = @animate for shift in shifts

    m = vec(getshift(Msub, shift))

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
           metricdisp=[:area, :shift,
                       :meanrate, :coh_at_zero], 
           scale_spg=2, 
           totalarea_aspect=1,
           metricfilter,
           srt=[:area,:bitsperspike],
           title_width=20,
           area_conditional_c)
Plot.save("fields_at_BEST_shift_srt=>bitsperspike")
