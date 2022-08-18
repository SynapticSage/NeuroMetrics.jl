"""
Correlation(τ) and Correlation(rateᵢⱼ)
"""

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
using Munge.spiking
using Filt
using DataFrames, DataFramesMeta
using Statistics, NaNStatistics
using Plots
using GLM
using StatsBase, StatsPlots
using ElectronDisplay
Plot.setfolder("timeshift","functional_connectivity")

@time spikes, beh, cells = Load.load("RY16", 36, data_source=["spikes","behavior", "cells"])
@time I = load_mains();
@time F = load_fields();
filt = get_filters()

datacut = :cue
spikes, beh = Load.load("RY16",36; data_source=["spikes","behavior"])
@time i = I[bestpartialmatch(keys(F), (;datacut, widths=5))];
@time f = F[bestpartialmatch(keys(F), (;datacut, widths=5))];
f = ShiftedFields(deepcopy(f))
unitshift = Timeshift.types.matrixform(f)
getshift(M::DimArray, s) = M[:, M.dims[2].==s];
shifts = collect(unitshift.dims[2])
push_celltable!( unitshift, cells, :unit, :area)
push_dims!(unitshift)
push_shiftmetric!(unitshift, best_tau!; metric=:bitsperspike)
#push_dims!(f, vec(getindex.(m, :coherence)); dim=:unit, metric=:coh_at_zero)

# Which cells pass our criteria?
metricfilter = x->
(x[:area] .=="CA1" .&& x[:meanrate] .> 0.01 && x[:meanrate] .< 7 .&& x[:maxrate] < 35) .||
(x[:area] .=="PFC" .&& x[:meanrate] .> 0.005 .&& x[:maxrate] < 100)
shift0 = filter(metricfilter, getshift(unitshift,0))
shift0 = shift0[sortperm(shift0[:unit])]
spikes = subset(spikes, :unit=>x->Utils.squeeze(any(x .∈ shift0[:unit]',dims=2)))
spikes = Utils.filtreg.filterAndRegister(beh,spikes; filters=filt[datacut], on="time",  transfer=["velVec"], filter_skipmissingcols=true)[2]

sort(unique(spikes.unit))

# ------------------------
# Firing rate correlations
# ------------------------

# Snage a rate matrix (time × neuron)
R = torate(spikes , beh; gaussian=0)
# Zscore
R = (R .- Array(mean(R, dims=:time)))./Array(std(R,dims=:time))

heatmap(R', clim=(-3, 3))
shift0_filtbyrate = shift0[vec(any(shift0[:unit] .∈ collect(R.dims[2])', dims=2))]

# Compute correlation
Rcor = cor(R, dims=:time)
heatmap(Array(Rcor), c=:vik, clim=(-1,1))

# ------------------------
# Unit τ differences Δ
# ------------------------

τ = NaN .* ones(size(Rcor))
(u1,u2,r) = first(zip(Rcor.dims[1], Rcor.dims[2], Rcor))
for (i1,i2) in Iterators.product(1:size(Rcor,1), 1:size(Rcor,2))
    u1, u2 = Rcor.dims[1].val[i1],
             Rcor.dims[2].val[i2]
    if u1 == u2
        continue
    end
    τ1, τ2 = shift0_filtbyrate[findfirst(shift0_filtbyrate[:unit] .== u1)][:bestshift_bitsperspike],
             shift0_filtbyrate[findfirst(shift0_filtbyrate[:unit] .== u2)][:bestshift_bitsperspike]
    τ[i1,i2] = τ1 - τ2
end


# ------------------------
# Rate corr ∝ Δτ
# ------------------------

using NaNStatistics
vecTau = vec(abs.(τ))
vecCor = vec(abs.(Rcor))
good = (!).(isnan.(vecTau))
vecTau, vecCor = vecTau[good], vecCor[good]
df = DataFrame(Dict(:vecTau=>vecTau, :vecCor=>vecCor))
df.rank = sortperm(df.vecCor)
df.bintau = Int.(floor.(Utils.norm_extrema(df.vecTau, [0,4])))
combine(groupby(df, :bintau),nrow,:vecTau=>mean)
#transform!(df, :vecCor => (x->log10.(x)) => :vecCor)

q = 95
@df df scatter(:bintau, :vecCor, label="samples", xlabel="δτ=τ₁-τ₂", ylabel="pearson")
dfs = combine(groupby(df, :bintau), :vecCor=>(x->quantile(x,q/100)))
@df dfs plot!(:bintau, :vecCor_function, label="Q_$q")
dfs = combine(groupby(df, :bintau), :vecCor=>mean)
@df dfs plot!(:bintau, :vecCor_mean, label="mean")
lm = fit(LinearModel, @formula(vecCor ~ vecTau + vecTau^2), df)
ŷ = predict(lm, df)
plot!(df.vecTau, ŷ, c=:black, title="Datacut=$datacut")
lm

q = 95
@df df scatter(:bintau, :rank, label="samples", xlabel="δτ=τ₁-τ₂", ylabel="pearson")
dfs = combine(groupby(df, :bintau), :rank=>(x->quantile(x,q/100)))
@df dfs plot!(:bintau, :rank_function, label="Q_$q")
dfs = combine(groupby(df, :bintau), :rank=>mean)
@df dfs plot!(:bintau, :rank_mean, label="mean")
lm = fit(LinearModel, @formula(rank ~ vecTau + vecTau^2), df)
ŷ = predict(lm, df)
plot!(df.vecTau, ŷ, c=:black, title="Datacut=$datacut")
lm

alpha=0.4
ecdfplot( @subset(df, :bintau .== 0).vecCor; fillrange=0, alpha, label="1st bin")
ecdfplot!(@subset(df, :bintau .== 1).vecCor; fillrange=0, alpha, label="2nd bin")
ecdfplot!(@subset(df, :bintau .== 2).vecCor; fillrange=0, alpha, label="3rd bin", legend=:outerbottomright)
ecdfplot!(@subset(df, :bintau .== 3).vecCor; alpha, fillrange=0, label="4th bin", legend=:outerbottomright)
ecdfplot!(@subset(df, :bintau .== 4).vecCor; alpha, fillrange=0, label="5th bin", legend=:outerbottomright)

@df df violin(:vecCor, group=:bintau)

import Gadfly
p = Gadfly.plot(x=vecTau, y=vecCor,
            Gadfly.Stat.binmean, Gadfly.Geom.point, Gadfly.Geom.line)


# ---------------
# XY CORRELATIONS
# ---------------

Plot.setfolder("timeshift", "xy-dist")
c = get.([ColorSchemes.vik],
        Utils.norm_extrema(shift0[:bestshift_bitsperspike], [0,1]))
task.plotboundary(taskdf, seriestype=:shape, alpha=0.2)
scatter!(eachcol(hcat(shift0[:centroid]...)')...; c, label="", framestyle=:none, title="Field centroids annotated betst τ")
Plot.save((;desc="Best_tau_of_field_centroids"))

   

# -------------------------
# FUTURE/PRESENT/PAST INFO
# -------------------------

# Create quintiles of bestshift and plot field averages
shift0[:bestshift_bins] =  Utils.binning.digitize(shift0[:bestshift_bitsperspike], 5)

averages_t0 = Dict(
    :standard=>[],
    :norm01 =>[],
    :normBSP=>[],
    :normB=>[],
    )
for bin = 1:5
    fields_in_bin = shift0[shift0[:bestshift_bins] .== bin]
    to_average = [f.rate for f in fields_in_bin]
    bitsperspike  = fields_in_bin[:bitsperspike]
    bits = fields_in_bin[:bitsperspike] .* fields_in_bin[:totalcount]
    push!(averages_t0[:standard], mean(to_average))
    push!(averages_t0[:norm01],   mean(map(x->Utils.nannorm_extrema(x, [0, 1]),
                                        to_average)))
    bitsperspike_norm(x::AbstractArray,bsp::Real) = Utils.nannorm_extrema(x, [0, bsp])
    push!(averages_t0[:normBSP],  mean(map(bitsperspike_norm,
                                        to_average, bitsperspike)))
    push!(averages_t0[:normB],  mean(map(bitsperspike_norm,
                                        to_average, bits)))
end

# Field averages, without scaling each neuron
plot([heatmap(aver, aspect_ratio=1, framestyle=:none) 
      for aver in averages_t0[:standard]]...)

# Field averages, scaling each neuron 0-1
plot([heatmap(aver, aspect_ratio=1, framestyle=:none) 
      for aver in averages_t0[:norm01]]...)
# Field averages, scaling each neuron by bits per spike
plot([heatmap(aver, aspect_ratio=1, framestyle=:none) 
      for aver in averages_t0[:normBSP]]...)

plot([heatmap(aver, aspect_ratio=1, framestyle=:none) 
      for aver in averages_t0[:normB]]...)


