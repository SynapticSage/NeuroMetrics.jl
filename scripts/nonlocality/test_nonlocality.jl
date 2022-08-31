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
using Munge.spiking
using Filt
using DataFrames, DataFramesMeta
using Statistics, NaNStatistics, StatsBase, StatsPlots, HypothesisTests, GLM
using Plots
#using ElectronDisplay
using LazySets
using Munge.timeshift: getshift
using Utils.stats: pfunc


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
datacut    = :all
datacutStr = string(datacut)

# Grab the key of interest given our datacut
@time i = I[bestpartialmatch(keys(F), (;datacut, widths=5, coactivity=nothing), nothing_means_removekey=true)];
@time f = F[bestpartialmatch(keys(F), (;datacut, widths=5, coactivity=nothing), nothing_means_removekey=true)];
f = ShiftedFields(deepcopy(f))
unitshift = Timeshift.types.matrixform(f)

# Setup a  shift-getting convenience method, the shifts, and a few metrics
shifts = collect(unitshift.dims[2])

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


function annotate_nonlocal_spikes!(spikes::DataFrame, unitshift::DimArray, shift::Union{<:Real, Symbol}=0;
                                    thresh=0.8, hullmax=1,
                                    n_examples=20, shuf=true)
    # Push details about fields and best shift
    push_celltable!(unitshift, cells, :unit, :area)
    push_dims!(unitshift)
    push_shiftmetric!(unitshift, best_tau!; metric=:bitsperspike)
    shift0 = unitshift[shift=At(shift)]
    # Add convex hull
    push_metric!(shift0, metrics.convexhull; thresh=0.80)
    hsg = [x[:convexhull][:hullseg_grid] for x in shift0];
    shift0[:hullseg_grid] = hsg;
    shift0 = DimArray(shift0, Dim{:unit}(shift0[:unit]));
    # Obtain nonlocal spikes
    P = []
    spikes.infield .= Vector{Union{Missing,Bool}}(undef, size(spikes,1))
    #trans(X) = [x[end:-1:begin] for x in X]
    for unit in shift0[:unit]
        inds = findall(spikes.unit .== unit .&& (!).(ismissing.(spikes.x)))
        if isempty(inds) || shift0[unit=At(unit)][:totalcount] == 0; 
        @info "skipping" unit
        continue; 
        else
        @info "processing" unit
        end
        spu = view(spikes, inds, :)
        histogram2d(spu.x, spu.y)
        cellspace = OrderedDict(r.time => Float32.(element(Singleton(r.x,r.y)))
                        for r in eachrow(spu))
        U = shift0[unit=At(unit)]
        hull = [VPolygon(U[:hullseg_grid][i])
                for i in 1:length(U[:hullseg_grid])][Utils.na, :]
        #hullT = [VPolygon(trans(U[:hullseg_grid][i]))
        #        for i in 1:length(U[:hullseg_grid])][Utils.na, :]
        if isempty(hull)
            @warn "no hull"
            continue
        else
            #@infiltrate
        end
        V = collect(values(cellspace))
        infield = V .∈ hull[:,1:hullmax]
        #infieldT = V .∈ hullT
        spu.infield = vec(any(infield,dims=2))
        title="unit=$unit,\nmean infield=$(round(mean(infield), sigdigits=2))"
        push!(P, Plots.plot(U; transpose=true))
        for h in 1:hullmax
        Plots.plot!(hull[h]; title)
        end
    end
    if shuf == true 
        Random.shuffle(P) 
    end
    Plots.plot(P[1:n_examples]...; 
    colorbar=false, titlefontsize=6, tickfontsize=6, framestyle=:none, size=(800,800))
end

annotate_nonlocal_spikes!(spikes, unitshift)
Plot.save((;desc="example of convex hulls"))

# Get likelihood of nonlocal
pfc_units = @subset(cells,:area.=="PFC").unit
R = Munge.spiking.torate(allspikes, beh)
pfc_units = intersect(pfc_units, collect(R.dims[2]))
infield = first(groupby(subset(spikes, :infield=>x->(!).(isnan.(x))) ,
                                        :infield))

# ------------
# Mean of rates, per cell
# ------------
infield_mean = DataFrame()
for infield in groupby(spikes, :infield)
    sample = [R[time=B, unit=At(pfc_units)]
                for B in Between.(infield.time.-0.02, infield.time.+0.02)]
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

ifnames =  OrderedDict(0.0 => :out_of_field, 1.0=>:infield)
transform!(infield_mean, DataFrames.All(), :infield => (x->[ifnames[xx] for xx in x]) => :infield_label)
@subset!(infield_mean, :infield .== 0 .|| :infield .== 1)
srt = SignedRankTest(@subset(infield_mean,:infield.==0).rate,
              @subset(infield_mean,:infield.==1).rate)
@df infield_mean boxplot(:infield, :rate, title="SignedRankPval=$(round(pvalue(srt),sigdigits=2)), with N=$(srt.n) PFC cells",
                         xticks=([0,1], collect(values(ifnames))))
@df infield_mean scatter!(:infield, :rate, group=:infield)

# ------------------------------
# Mean of rates, per cell x time
# ------------------------------
infield_unit = DataFrame()
for infield in groupby(spikes, :infield)
    sample = [R[time=B, unit=At(pfc_units)]
                for B in Between.(infield.time.-0.02, infield.time.+0.02)]
    sample = [s for s in sample if !isempty(s)]
    @time sample = vcat(sample...)
    unit = repeat(vec(pfc_units)', size(sample,1))
    append!(infield_unit,
    DataFrame(OrderedDict(
        :unit => vec(unit),
        :infield => infield.infield[1],
        :rate => vec(sample) ))
        )
    #push!(infield_median, infield.infield[1] => rate_median)
end

@subset!(infield_unit, :infield .== 0 .|| :infield .== 1)
transform!(infield_unit, DataFrames.All(), :infield => (x->[ifnames[xx] for xx in x]) => :infield_label)
srt_units

srt = OneWayANOVATest(@subset(infield_unit,:infield.==0).rate,
              @subset(infield_unit,:infield.==1).rate)
@df infield_unit boxplot(:infield, :rate, title="TTest=$(round(pvalue(srt),sigdigits=2)) difference of individual cells",
                         xticks=([0,1], collect(values(ifnames))), outliers=false)
Plot.save((;desc="difference of individual cells"))


srt   = OrderedDict()
diffs = OrderedDict()
@showprogress for U in groupby(infield_unit, :unit)
    push!(srt,
            U.unit[1] => UnequalVarianceTTest(@subset(U,:infield.==0).rate, @subset(U,:infield.==1).rate))
    push!(diffs,
          U.unit[1] => mean(@subset(U,:infield.==0).rate) - mean(@subset(U,:infield.==1).rate))

end

pv_higherout = pvalue.(values(srt)) .< 0.05/length(srt) .&& values(diffs) .> 0
pv_higherin  = pvalue.(values(srt)) .< 0.05/length(srt) .&& values(diffs) .< 0
pv_nonsig    = pvalue.(values(srt)) .> 0.05/length(srt)
res =[collect(values(diffs)) pvalue.(values(srt))]
bar([0, 1, 2], mean.([pv_nonsig, pv_higherin, pv_higherout]),
    xticks=([0, 1, 2], ["not bonferroni\nsig", "sig in-field", "sig out-field"]),
    ylabel="Percent", title="Nonlocality: PFC firing conditioned on CA1 field")
Plot.save((;desc="bonferroni sig cells"))

# -----------------------------------------------------
# Split by cuemem
# -----------------------------------------------------
Utils.filtreg.register(beh,spikes, on="time", transfer=["cuemem"])

# RATE MOD
infield_cuememunit = DataFrame()
@showprogress for group in groupby(dropmissing(spikes,:cuemem), [:infield,:cuemem])
    @info "samples" size(group.time)
    sample = [R[time=B, unit=At(pfc_units)]
                for B in Between.(group.time.-0.02, group.time.+0.02)]
    sample = [s for s in sample if !isempty(s)]
    @time sample = vcat(sample...)
    unit = repeat(vec(pfc_units)', size(sample,1))
    @infiltrate
    append!(infield_cuememunit,
    DataFrame(OrderedDict(
        :unit => vec(unit),
        :infield => group.infield[1],
        :cuemem => group.cuemem[1],
        :rate => vec(sample) ))
        )
end

clab = OrderedDict(-1 => "nontask", 0 => "cue", 1=> "mem")

D = DataFrame()
@time @showprogress for group in groupby(infield_cuememunit,[:cuemem])
     outrate= @subset(group,:infield.==0).rate
     inrate= @subset(group,:infield.==1).rate
     test =  UnequalVarianceTTest(outrate, inrate)
     pval = pvalue(test)
     append!(D, OrderedDict(
        :cuemem => group.cuemem[1],
        :cuemem_label => clab[group.cuemem[1]],
        :diff => mean(outrate) - mean(inrate),
        :pval => pval,
        :test => test,
       ))
end
transform!(D, DataFrames.All(), :pval => pfunc => :pval_label)
@df D bar(:cuemem, :diff, xticks=([-1,0,1], collect(values(clab))) )
markers = text.(D.pval_label, :white)
@df D annotate!(:cuemem, :diff * 1.025, markers)

Du = DataFrame()
@time @showprogress for group in groupby(infield_cuememunit,[:cuemem,:unit])
     outrate= @subset(group,:infield.==0).rate
     inrate= @subset(group,:infield.==1).rate
     test =  UnequalVarianceTTest(outrate, inrate)
     pval = pvalue(test)
     append!(D, OrderedDict(
        :cuemem => group.cuemem[1],
        :cuemem_label => clab[group.cuemem[1]],
        :unit => group.unit[1],
        :diff => mean(outrate) - mean(inrate),
        :pval => pval,
        :test => test,
       ))
end

transform!(D, DataFrames.All(), :pval => pfunc => :pval_label)
@df D bar(:cuemem, :diff, xticks=([-1,0,1], collect(values(clab))) )
markers = text.(D.pval_label, :white)
@df D annotate!(:cuemem, :diff * 1.025, markers)



# SHEER NUMBER OF OUT OF FIELDERS
if_counts = combine(groupby(spikes, [:cuemem, :infield]), nrow)
subset!(if_counts, :infield => x->(!).(isnan.(x)))
if_perc = combine(groupby(if_counts, [:cuemem]), :infield, :nrow => (x->x./(sum(x))) => :perc)
iflab = Dict(0 => "out", 1=> "in")
clab = Dict(-1 => "nontask", 0 => "cue", 1=> "mem")
transform!(if_perc, DataFrames.All(), [:cuemem, :infield] => ((x,y)->(getindex.([clab], x) .* " " .* getindex.([iflab], y))) .=> :cuemem_inoutfield)
sort!(if_perc, :infield)
@df if_perc bar(:cuemem_inoutfield, :perc, group=:infield)

Plot.save((; desc="sheer fraction of outfielder spikes by cue and mem and nontask"))


# -----------------------------------------------------
# Reverse the direction: CA1 out of field per PFC spike
# -----------------------------------------------------


Plot.setfolder("nonlocality", "pfc_rates")
Rpfc =  Munge.spiking.torate(@subset(allspikes, :area .== "PFC"), beh)
Dpfc = DataFrame([0;Rpfc], :auto)
Dpfc.cuemem = beh.cuemem
pfc_cell_means = Matrix(combine(groupby(Dpfc, "cuemem"), 1:21 .=> mean))[:,2:end]
pfc_cell_means./maximum(pfc_cell_means,dims=1)
heatmap(["Nontask","Cue","Memory"], 1:size(pfc_cell_means,2), pfc_cell_means./maximum(pfc_cell_means,dims=1))
Plot.save("Average PFC cell firing rates per task type")

bar(["Nontask","Cue","Memory"], mean(pfc_cell_means./maximum(pfc_cell_means,dims=1), dims=2))
Plot.save("Average of PFC cell firing, mean(norm0)")
bar(["Nontask","Cue","Memory"], mean(pfc_cell_means, dims=2))
Plot.save("Average of PFC cell firing, mean(rate)")

