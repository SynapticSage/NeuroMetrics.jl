"""
Correlation(τ) and Correlation(rateᵢⱼ)

TODOs
----
- Higher resolution for XY analyses

"""

using DrWatson
using Timeshift
using Timeshift.types
using Timeshift.shiftmetrics
using Field.metrics
using Plot
using Plot.receptivefield
using Utils.namedtup
using DataStructures: OrderedDict
using DimensionalData
import DimensionalData: Between
using ProgressMeter
import Plot
using Munge.spiking
using Filt
using DataFrames, DataFramesMeta
using Statistics, NaNStatistics
using Plots
using GLM
using StatsBase, StatsPlots
#using ElectronDisplay
using ProgressMeter
using LazySets
Plot.setfolder("nonlocality")

# Get our dataframes of interest and the field objects
@time I = load_mains();
@time F = load_fields();
filt = get_filters()

@time spikes, beh, cells = Load.load("RY16", 36, data_source=["spikes","behavior", "cells"])
beh, spikes = Utils.filtreg.register(beh,spikes, on="time", transfer=["x","y"])
allspikes = copy(spikes)
spikes.infield = fill(NaN, size(spikes.time))

# Setup our cell filters
metricfilters = Dict(
    :CA1=> x->
    x[:area] .=="CA1" .&& x[:meanrate] .> 0.01 && x[:meanrate] .< 7 .&& x[:maxrate] < 35 .&& x[:coherence] > 0.5,
    :PFC => x->
    x[:area] .=="PFC" .&& x[:meanrate] .> 0.005 .&& x[:maxrate] < 100,
   )
metricfilters[:CA1PFC] = x-> metricfilters[:CA1](x) || metricfilters[:PFC](x)


# Which cut of data for fields and spikes?
datacut = :all
datacutStr = string(datacut)

# Grab the key of interest given our datacut
@time i = I[bestpartialmatch(keys(F), (;datacut, widths=5))];
@time f = F[bestpartialmatch(keys(F), (;datacut, widths=5))];
f = ShiftedFields(deepcopy(f))
unitshift = Timeshift.types.matrixform(f)

# Setup a  shift-getting convenience method, the shifts, and a few metrics
getshift(arrayOfFields::DimArray, s) = arrayOfFields[:, arrayOfFields.dims[2].==s];
shifts = collect(unitshift.dims[2])
push_celltable!( unitshift, cells, :unit, :area)
push_dims!(unitshift)
push_shiftmetric!(unitshift, best_tau!; metric=:bitsperspike)

# Which cells pass our criteria?
region = :CA1
metricfilter = metricfilters[region]

# Get filtered shift=0 fields
shift0 = filter(metricfilter, getshift(unitshift,0))
shift0 = shift0[sortperm(shift0[:unit])]

# Subset spikes by our filtration schema
spikes = subset(spikes, :unit=>x->Utils.squeeze(any(x .∈ shift0[:unit]',dims=2)))
spikes = Utils.filtreg.filterAndRegister(beh,spikes; filters=filt[datacut], on="time",  transfer=["velVec"], filter_skipmissingcols=true)[2]
sort(unique(spikes.unit))

# Add convex hull
push_metric!(shift0, metrics.convexhull; thresh=0.85)
hsg = [x[:convexhull][:hullseg_grid] for x in shift0];
shift0[:hullseg_grid] = hsg;
shift0 = DimArray(shift0, Dim{:unit}(shift0[:unit]));

# Obtain nonlocal spikes
P = []
for unit in shift0[:unit]
    @info "unit" unit
    inds = findall(spikes.unit .== unit .&& (!).(ismissing.(spikes.x)))
    if isempty(inds); continue; end
    spu = view(spikes, inds, :)
    histogram2d(spu.x, spu.y)
    cellspace = OrderedDict(r.time => Float32.(element(Singleton(r.x,r.y)))
                    for r in eachrow(spu))
    U = shift0[unit=At(unit)]
    hull = VPolygon(U[:hullseg_grid][1])
    V = collect(values(cellspace))
    infield = V .∈ [hull]
    spu.infield = infield
    title="unit=$unit,\nmean infield=$(round(mean(infield), sigdigits=2))"
    push!(P, plot(U; transpose=true))
    plot!(hull; title)
end
plot(P[1:20]...; colorbar=false, titlefontsize=6, tickfontsize=6, framestyle=:none)


# Get likelihood of nonlocal
pfc_units = @subset(cells,:area.=="PFC").unit
R = Munge.spiking.torate(allspikes, beh)
pfc_units = intersect(pfc_units, collect(R.dims[2]))
infield = first(groupby(subset(spikes, :infield=>x->(!).(isnan.(x))) ,
                                        :infield))

# Mean rates
infield_mean = DataFrame()
for infield in groupby(spikes, :infield)
    sample = [R[time=B, unit=At(pfc_units)] for B in Between.(infield.time.-0.02, infield.time.+0.02)]
    sample = [s for s in sample if !isempty(s)]
    @time sample = vcat(sample...)
    rate_mean   = mean(sample, dims=1)
    #rate_median = median(sample[sample.!=0])
    append!(infield_mean,
    DataFrame(OrderedDict(
    :unit => vec(pfc_units),
    :infield => infield.infield[1],
    :rate => vec(rate_mean))))
    #push!(infield_median, infield.infield[1] => rate_median)
end
@subset!(infield_mean, :infield .== 0 .|| :infield .== 1)
#@df subset(infield_mean, :infield => x->x .== 0) boxplot(:infield, :rate)
@df infield_mean boxplot(:infield, :rate)

