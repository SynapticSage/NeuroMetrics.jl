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
using Statistics, NaNStatistics, StatsBase, StatsPlots, HypothesisTests, GLM
using Plots
#using ElectronDisplay
using LazySets
using Munge.timeshift: getshift
using Utils.stats: pfunc
clab = OrderedDict(-1 => "nontask", 0 => "cue", 1=> "mem", missing=>"sleep")
isonames =  OrderedDict(false => :adjacent, true=>:isolated)
filt_desc = OrderedDict(:all => "> 2cm/s")
save_kws = (;pfc_rate_analy=true)


# Acquire data
@time spikes, beh, cells = Load.load("RY16", 36, data_source=["spikes","behavior", "cells"])
beh, spikes = Utils.filtreg.filterAndRegister(beh, spikes, on="time", transfer=["x","y"], filters=filt[:all])
allspikes = copy(spikes)
beh2 = Load.load_behavior("RY16",36)

# Acquire LFP and isolated spikes
lfp = Load.load_lfp("RY16", 36, tet=5);
lfp.time = lfp.time .- Load.min_time_records[1]
lfp = Munge.lfp.annotate_cycles(lfp)
sp = @subset(spikes, :tetrode .== 6);
Munge.spiking.isolated(sp, lfp)

# Split by isolated spikes and discover the fraction of isolated spikes
iso_sum = combine(groupby(spikes, [:area, :cuemem]), :isolated => mean, (x->nrow(x)))

# Calculate time animal spends in each cuemem segment
task_pers = Table.get_periods(beh2, [:traj, :cuemem], timefract= :velVec => x->abs(x) > 2)

# Total that time and register that column to the isolation summary
task_pers = combine(groupby(task_pers, [:cuemem]), [:δ,:frac] => ((x,y)->sum(x.*y)) => :timespent)
Utils.filtreg.register(task_pers, iso_sum, on="cuemem", transfer=["timespent"])

# Acqruire events per time as events  / time spent
iso_sum.events_per_time = iso_sum.x1 ./ (iso_sum.timespent)
iso_sum.cuearea = iso_sum.area .* "\n" .* getindex.([clab], iso_sum.cuemem)
iso_sum.cmlab = getindex.([clab], iso_sum.cuemem)
iso_sum.isolated_events_per_time = iso_sum.isolated_mean .* iso_sum.events_per_time

ord = Dict("nontask"=>1,"cue"=>2,"mem"=>3)
iso_sum = sort(iso_sum, [DataFrames.order(:cmlab, by=x->ord[x]),:cuearea])

Plot.setfolder("nonlocality","MUA and isolated MUA")
# =========PLOTS =======================================
kws=(;legend_position=Symbol("outerbottomright"))
@df @subset(iso_sum,:area.=="CA1") bar(:cmlab, :timespent, ylabel="Time spent", group=:cuemem; kws...)
Plot.save((;desc="time spent"))
@df iso_sum bar(:cuearea, :events_per_time, ylabel="MUA events per second\n$(filt_desc[:all])", group=:cuemem; kws...)
Plot.save((;desc="MUA per second"))
@df iso_sum bar(:cuearea, :isolated_mean, ylabel="Isolated Spikes (sign of CA1-PFC interaction)\n$(filt_desc[:all])", group=:cmlab; kws...)
Plot.save((;desc="fraction of isolated spikes"))
@df iso_sum bar(:cuearea, :isolated_events_per_time, ylabel="Isolated MUA × sec⁻1\n(sign of CA1-PFC interaction)\n$(filt_desc[:all])", group=:cuemem; kws...)
Plot.save((;desc="isolated spikes per second"))
# ======================================================

Plot.setfolder("nonlocality", "isolated_ca1_rate_pfc")


# Isolated spikees by ISI
# -----------------------
Munge.spiking.nextandprev!(spikes)
perc = round(mean(spikes.neard .> .4)*100,sigdigits=2)
histogram(log10.(spikes.neard),label="nearest spike")
vline!([log10(.400)], label="thresh", title="$perc % isolated")
Plot.save("nearest spike isolation")


pfc_units = @subset(cells,:area.=="PFC").unit
R = Munge.spiking.torate(allspikes, beh)
pfc_units = intersect(pfc_units, collect(R.dims[2]))

isolated = last(groupby(subset(spikes, :isolated=>x->(!).(isnan.(x))) ,
                                        :isolated))
@assert all(isolated.isolated .== true)

# ----------------------
"Mean of rates, per cell"
# ----------------------

isolated_mean = DataFrame()
for isolated in groupby(spikes, :isolated)
    sample = [R[time=B, unit=At(pfc_units)]
                for B in Between.(isolated.time.-0.02, isolated.time.+0.02)]
    sample = [s for s in sample if !isempty(s)]
    @time sample = vcat(sample...)
    rate_mean   = mean(sample, dims=1)
    #rate_mean = mean(sample[sample.!=0])
    append!(isolated_mean,
    DataFrame(OrderedDict(
    :unit => vec(pfc_units),
    :isolated => isolated.isolated[1],
    :rate => vec(rate_mean))))
    #push!(isolated_mean, isolated.isolated[1] => rate_mean)
end

transform!(isolated_mean, DataFrames.All(), :isolated => (x->[isonames[xx] for xx in x]) => :isolated_label)
@subset!(isolated_mean, :isolated .== 0 .|| :isolated .== 1)
srt = SignedRankTest(@subset(isolated_mean,:isolated.==0).rate,
              @subset(isolated_mean,:isolated.==1).rate)
@df isolated_mean boxplot(:isolated, :rate, title="SignedRankPval=$(round(pvalue(srt),sigdigits=2)), with N=$(srt.n) PFC cells",
                         xticks=([0,1], collect(values(isonames))))
@df isolated_mean scatter!(:isolated, :rate, group=:isolated)
Plot.save((;save_kws...,desc="meanmean_pfc_cell_rate"))

# ------------------------------
# Mean of rates, per cell x time
# ------------------------------
isolated_unit = DataFrame()
for cell_isolation_group in groupby(spikes, :isolated)
    sample = [R[time=B, unit=At(pfc_units)]
                for B in Between.(cell_isolation_group.time.-0.02,
                                  cell_isolation_group.time.+0.02)]
    sample = [s for s in sample if !isempty(s)]
    @time sample = vcat(sample...)
    unit = repeat(vec(pfc_units)', size(sample,1))
    append!(isolated_unit,
    DataFrame(OrderedDict(
        :unit     => vec(unit),
        :isolated => cell_isolation_group.isolated[1],
        :rate     => vec(sample) ))
        )
end

@subset!(isolated_unit, :isolated .== 0 .|| :isolated .== 1)
transform!(isolated_unit, DataFrames.All(), :isolated => (x->[isonames[xx] for xx in x]) => :isolated_label)
#srt_units

srt = UnequalVarianceTTest(@subset(isolated_unit,:isolated.==0).rate,
              @subset(isolated_unit,:isolated.==1).rate)
@df isolated_unit boxplot(:isolated, :rate, title="Firing rate samples per spike\nTTest=$(round(pvalue(srt),sigdigits=2)) difference of individual cells",
                         xticks=([0,1], collect(values(isonames))), outliers=false)
Plot.save((;save_kws...,desc="firing rate samples per spike"))

"""
DO SAME split by cue and memory!

(STOPPED UPDATING HERE)
"""
#TODO

srt   = OrderedDict()
diffs = OrderedDict()
@showprogress for U in groupby(isolated_unit, :unit)
    push!(srt,
            U.unit[1] => UnequalVarianceTTest(@subset(U,:isolated.==0).rate,
                                              @subset(U,:isolated.==1).rate))
    push!(diffs,
          U.unit[1] => mean(@subset(U,:isolated.==0).rate) - mean(@subset(U,:isolated.==1).rate))

end

pv_higheradjacent = pvalue.(values(srt)) .< 0.05/length(srt) .&& values(diffs) .> 0
pv_higherisolated  = pvalue.(values(srt)) .< 0.05/length(srt) .&& values(diffs) .< 0
pv_nonsig    = pvalue.(values(srt)) .> 0.05/length(srt)
res =[collect(values(diffs)) pvalue.(values(srt))]
bar([0, 1, 2], mean.([pv_nonsig, pv_higheradjacent, pv_higherisolated]),
    xticks=([0, 1, 2], ["not bonferroni\nsig", "cells sig adjacent", "cells sig isolated"]),
    ylabel="Percent", title="Nonlocality: PFC firing conditioned on CA1 field")

Plot.save((;save_kws...,desc="bonferroni sig per cell, firing rate samples per spike"))

# -----------------------------------------------------
# Split by cuemem
# -----------------------------------------------------
Utils.filtreg.register(beh,spikes, on="time", transfer=["cuemem"])

# RATE MOD
isolated_cuememunit = DataFrame()
@showprogress for group in groupby(dropmissing(spikes,:cuemem), [:isolated,:cuemem])
    @info "samples" size(group.time)
    sample = [R[time=B, unit=At(pfc_units)]
                for B in Between.(group.time.-0.02, group.time.+0.02)]
    sample = [s for s in sample if !isempty(s)]
    @time sample = vcat(sample...)
    unit = repeat(vec(pfc_units)', size(sample,1))
    append!(isolated_cuememunit,
    DataFrame(OrderedDict(
        :unit => vec(unit),
        :isolated => group.isolated[1],
        :cuemem => group.cuemem[1],
        :rate => vec(sample) ))
        )
end

clab = OrderedDict(-1 => "nontask", 0 => "cue", 1=> "mem")

D = DataFrame()
@time @showprogress for group in groupby(isolated_cuememunit,[:cuemem])
     adjacent_spikes = @subset(group,:isolated.==0).rate
     iso_spikes      = @subset(group,:isolated.==1).rate
     test =  UnequalVarianceTTest(adjacent_spikes, iso_spikes)
     pval = pvalue(test)
     append!(D, OrderedDict(
        :cuemem => group.cuemem[1],
        :cuemem_label => clab[group.cuemem[1]],
        :diff            => mean(iso_spikes) - mean(adjacent_spikes),
        :iso_spikes      => mean(iso_spikes),
        :adjacent_spikes => mean(adjacent_spikes),
        :pval => pval,
        :test => test,
       ))
end
Du = DataFrame()
@time @showprogress for group in groupby(isolated_cuememunit,[:cuemem,:unit])
     adjacent_spikes = @subset(group,:isolated.==0).rate
     iso_spikes      = @subset(group,:isolated.==1).rate
     test =  UnequalVarianceTTest(adjacent_spikes, iso_spikes)
     pval = pvalue(test)
     append!(Du, OrderedDict(
        :cuemem => group.cuemem[1],
        :cuemem_label => clab[group.cuemem[1]],
        :unit => group.unit[1],
        :diff => mean(iso_spikes) - mean(adjacent_spikes),
        :iso_spikes      => mean(iso_spikes),
        :adjacent_spikes => mean(adjacent_spikes),
        :pval => pval,
        :test => test,
       ))
end

D.clabel = getindex.([clab], D.cuemem)

transform!(D, Colon(), :pval => Utils.pfunc => :pval_label)
#xticks=([-1,0,1], collect(values(clab)
@df D bar(:clabel, :diff, group=:clabel)
markers = text.(D.pval_label, :white)
@df D annotate!(:clabel, :diff * 1.025, markers)
#transform!(Du, DataFrames.All(), :pval => Utils.pfunc => :pval_label)
Plot.save((;save_kws..., desc="diff of iso adjacent in each state"))

@df Du scatter!(:cuemem .+ 0.05 .* randn(size(:cuemem)) , :diff, xticks=([-1,0,1], collect(values(clab))) )

#markers = text.(D.pval_label, :white)
#@df Du annotate!(:cuemem, :diff * 1.025, markers)

# SHEER NUMBER OF OUT OF FIELDERS
if_counts = combine(groupby(spikes, [:cuemem, :isolated]), nrow)
subset!(if_counts, :isolated => x->(!).(isnan.(x)))
if_perc = combine(groupby(if_counts, [:cuemem]), :isolated, :nrow => (x->x./(sum(x))) => :perc)
iflab = Dict(0 => "out", 1=> "in")
clab = Dict(-1 => "nontask", 0 => "cue", 1=> "mem")
transform!(if_perc, DataFrames.All(), [:cuemem, :isolated] => ((x,y)->(getindex.([clab], x) .* " " .* getindex.([iflab], y))) .=> :cuemem_inoutfield)
sort!(if_perc, :isolated)
@df subset(if_perc, :isolated => (x->x.!=false)) bar(:cuemem_inoutfield, :perc, group=:isolated)

Plot.save((;save_kws..., desc="cellcounts sig"))


# -----------------------------------------------------
# PFC firign rates norm by 0-1
# -----------------------------------------------------
Plot.setfolder("nonlocality", "isolated_ca1_pfc_rates")
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

# -----------------------------------------------------
# PFC firing rates norm by zscore
# -----------------------------------------------------
Plot.setfolder("nonlocality", "isolated_ca1_pfc_rates")
Rpfc =  Munge.spiking.torate(@subset(allspikes, :area .== "PFC"), beh)
Rpfc =  (Rpfc.-Matrix(mean(Rpfc,dims=1)))./Matrix(std(Rpfc,dims=1))
Dpfc = DataFrame([0;Rpfc], :auto)
Dpfc.cuemem = beh.cuemem
pfc_cell_means = Matrix(combine(groupby(Dpfc, "cuemem"), 1:21 .=> mean))[:,2:end]
pfc_cell_means./maximum(pfc_cell_means,dims=1)
heatmap(["Nontask","Cue","Memory"], 1:size(pfc_cell_means,2), pfc_cell_means./maximum(pfc_cell_means,dims=1))
Plot.save("Average PFC cell firing rates per task type")

bar(["Nontask","Cue","Memory"], mean(pfc_cell_means, dims=2))

Plot.save("Average of PFC cell firing, mean(norm0)")
bar(["Nontask","Cue","Memory"], mean(pfc_cell_means, dims=2))
Plot.save("Average of PFC cell firing, mean(rate)")

