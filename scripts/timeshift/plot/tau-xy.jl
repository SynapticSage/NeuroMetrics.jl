
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

task = Load.load_task("RY16",36)

# Setup basic metrics, like the area and dimensions of the dimarray
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
Msub0 = getshift(Msub, 0)
Msub = filter(metricfilter, Msub)
Plot.setfolder("timeshift","population","tau-xy")


good = findall((!).(isnan.(Msub0[:bestshift_bitsperspike])))
m = Msub0[good]
c = get.([ColorSchemes.vik], Utils.norm_extrema(m[:bestshift_bitsperspike], [0,1]))
scatter(eachcol(hcat(m[:centroid]...)')...)
scatter(eachcol(hcat(Msub0[:centroid]...)')...; c)
plotboundary!(task)



