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
filt = Filt.get_filters()
datacut = :all

# ===================
# ACQUIRE DATA
# ===================
# Acquire data
@time spikes, beh, cells = Load.load("RY16", 36, data_source=["spikes","behavior", "cells"])
beh, spikes = Utils.filtreg.filterAndRegister(beh, spikes, on="time", transfer=["x","y","cuemem"], filters=filt[datacut], filter_skipmissingcols=true)
allspikes = copy(spikes)
beh2 = Load.load_behavior("RY16",36)

# Acquire LFP and isolated spikes
lfp = Load.load_lfp("RY16", 36, tet=5);
lfp.time = lfp.time .- Load.min_time_records[1]
lfp = Munge.lfp.annotate_cycles(lfp)
#sp = @subset(spikes, :tetrode .== 6);

F = load_fields()
@time f = F[bestpartialmatch(keys(F), (;datacut, widths=5,coactivity=nothing), nothing_means_removekey=true)];
f = ShiftedFields(deepcopy(f))
unitshift = Timeshift.types.matrixform(f)
annotate_nonlocal_spikes!(spikes, unitshift, 0)

# ===================
# OUT OF FIELD SPIKES
# ===================
# Setup a  shift-getting convenience method, the shifts, and a few metrics
shifts = collect(unitshift.dims[2])
push_celltable!( unitshift, cells, :unit, :area)
push_dims!(unitshift)
push_shiftmetric!(unitshift, best_tau!; metric=:bitsperspike)

# Which cells pass our criteria?
#region = :CA1
#metricfilter = metricfilters[region]
#
## Get filtered shift=0 fields
#shift0 = filter(metricfilter, getshift(unitshift,0))
#shift0 = shift0[sortperm(shift0[:unit])]

# Subset spikes by our filtration schema
#spikes = subset(spikes, :unit=>x->Utils.squeeze(any(x .∈ shift0[:unit]',dims=2)))
#spikes = Utils.filtreg.filterAndRegister(beh,spikes; filters=filt[datacut], on="time",  transfer=["velVec"], filter_skipmissingcols=true)[2]
#sort(unique(spikes.unit))


# ===================
# ISOLATED SPIKING
# ===================
#Munge.spiking.isolated(sp, lfp, include_samples=true)
Munge.spiking.isolated(spikes, lfp, include_samples=false)


# Split by isolated spikes and discover the fraction of isolated spikes
function get_isolation_summary(spikes,split=[:cuemem])
    iso_sum = combine(groupby(dropmissing(spikes,:isolated), [:area, split...]), [:isolated,:nearestcyc,:meancyc] .=> mean, (x->nrow(x)))
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
end

"""
Title: Sheer diffferent rate of events (adjacent/isolated) over cue, mem, nontask
"""
iso_sum = get_isolation_summary(spikes)
sort!(iso_sum, [:area, :cuemem])

pfc_units = @subset(cells,:area.=="PFC").unit
R = Munge.spiking.torate(allspikes, beh)
pfc_units = intersect(pfc_units, collect(R.dims[2]))

isolated = last(groupby(subset(spikes, :isolated=>x->(!).(isnan.(x))) ,
                                        :isolated))
@assert all(isolated.isolated .== true)

# =========PLOTS =======================================
Plot.setfolder("nonlocality","MUA and isolated MUA")
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


Load.register(cells, spikes, on="unit", transfer=["meanrate"])
spikes.interneuron = spikes.meanrate .> 6
iso_sum_celltype = get_isolation_summary(spikes,[:cuemem, :interneuron])
sort!(iso_sum_celltype, [:area, :interneuron, :cuemem])

# =========PLOTS =======================================
@df iso_sum_celltype bar(:cuearea, :events_per_time, ylabel="MUA events per second\n$(filt_desc[:all])", group=:interneuron; kws...)
Plot.save((;desc="MUA per second, celltype"))
@df iso_sum_celltype bar(:cuearea, :isolated_mean, ylabel="Isolated Spikes (sign of CA1-PFC interaction)\n$(filt_desc[:all])", group=:interneuron; kws...)
Plot.save((;desc="fraction of isolated spikes, celltype"))
@df iso_sum_celltype bar(:cuearea, :isolated_events_per_time, ylabel="Isolated MUA × sec⁻1\n(sign of CA1-PFC interaction)\n$(filt_desc[:all])", group=:interneuron; kws...)
Plot.save((;desc="isolated spikes per second, celltype"))
# ======================================================

"""
Period-wise calculations
"""

Load.register(beh, spikes, on="time", transfer=["period","correct"])
iso_sum_celltype_per = get_isolation_summary(spikes, [:cuemem, :interneuron, :period, :correct])
sort!(iso_sum_celltype_per, [:area, :interneuron, :cuemem, :period])
@subset!(iso_sum_celltype_per, (:cuemem .== -1 .&& :correct .== -1) .||
                                (:cuemem .== 0 .&& :correct .!= -1) .||
                                (:cuemem .== 1 .&& :correct .!= -1))


# =========PLOTS =======================================
# -------- mua per second ------------------------------
Plot.setfolder("nonlocality","MUA and isolated MUA", "period-wise")
@df @subset(iso_sum_celltype_per, :interneuron .== false) scatter(:cuemem .+ (randn(size(:cuemem)) .* 0.1), :events_per_time, group=:correct, xtick=([-1,0,1],["nontask","cue","mem"]),title="MUA events/time", xlabel="task", ylabel="events × time\$^{-1}\$, pyr cells", legend_title="correct/error", legend_position=:outerbottomright)
Plot.save((;desc="mua-per-time",group=:correct,pyr=true))

@df @subset(iso_sum_celltype_per, :interneuron .== false) scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :events_per_time, group=:correct, alpha=0.5, xtick=([-1,0,1],["nontask","cue","mem"]),title="MUA events/time", xlabel="task", ylabel="events × time\$^{-1}\$, pyr cells", legend_title="correct/error", legend_position=:outerbottomright)
Plot.save((;desc="mua-per-time",group=:correct,pyr=true,sep_groups=true))

@df iso_sum_celltype_per scatter(:cuemem .+ (randn(size(:cuemem)) .* 0.1), :events_per_time, group=:interneuron, alpha=0.5, xtick=([-1,0,1],["nontask","cue","mem"]),title="MUA events/time", xlabel="task", ylabel="events × time\$^{-1}\$", legend_title="pyr/int", legend_position=:outerbottomright)
Plot.save((;desc="mua-per-time",group=:interneuron))

@df iso_sum_celltype_per scatter(:cuemem .+ 0.25.*(:interneuron .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :events_per_time, group=:interneuron, alpha=0.5, xtick=([-1,0,1],["nontask","cue","mem"]),title="MUA events/time", xlabel="task", ylabel="events × time\$^{-1}\$", legend_title="pyr/int", legend_position=:outerbottomright)
Plot.save((;desc="mua-per-time",group=:interneuron,group_sep=true))

# -------- percent isolated spikes ---------------------
@df @subset(iso_sum_celltype_per, :interneuron .== false) scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :isolated_mean, group=:correct, alpha=0.5, xtick=([-1,0,1],["nontask","cue","mem"]),title="Isolation", xlabel="task", ylabel="%percent isolated spikes", legend_title="correct/error", legend_position=:outerbottomright)
Plot.save((;desc="isolation", group=:correct, group_sep=true, ))
@df iso_sum_celltype_per scatter(:cuemem .+ 0.25.*(:interneuron .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :isolated_mean, group=:interneuron, alpha=0.5)
Plot.save((;desc="isolation", group=:interneuron, group_sep=true))

# -------- isolated per second spikes ------------------
@df @subset(iso_sum_celltype_per, :interneuron .== false) scatter(:cuemem .+ (randn(size(:cuemem)) .* 0.1), :isolated_events_per_time, group=:correct)
Plot.save((;desc="isolated-MUA-per-time", group=:correct, group_sep=true))
@df @subset(iso_sum_celltype_per, :interneuron .== false) scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :isolated_events_per_time, group=:correct, alpha=0.5)
Plot.save((;desc="isolated-MUA-per-time", group=:interneuron, group_sep=true))

# -------- multi - metric : one predict other ? --------
anim = @animate for i in 1:360
    p1=@df @subset(iso_sum_celltype_per, :interneuron .== false, :cuemem .== 0) scatter(:isolated_mean, :isolated_events_per_time, :events_per_time,
                group=:correct, alpha=0.5, camera=(i,30), 
                title="pyr cue",
                xlabel="fraction(isolated)",ylabel="isolated multiunit/s", zlabel="multiunit/s")
    p2=@df @subset(iso_sum_celltype_per, :interneuron .== false, :cuemem .== 1) scatter(:isolated_mean, :isolated_events_per_time, :events_per_time,
                group=:correct, alpha=0.5, camera=(i,30),
                title="pyr mem",
                xlabel="fraction(isolated)", ylabel="isolated multiunit/s", zlabel="multiunit/s")
    plot(p1,p2)
end
gif(anim, plotsdir(Plot.folder_args..., "separatrix-three-vars-gif_pyr_facet=cuemem_group=correct.gif"); loop=1)

anim = @animate for i in 1:360
    p1=@df @subset(iso_sum_celltype_per, :interneuron .== false, :cuemem .== 1) scatter(:isolated_mean, :isolated_events_per_time, :events_per_time, group=:correct, alpha=0.5, camera=(i,30), xlabel="fraction(isolated)",ylabel="isolated multiunit/s", zlabel="multiunit/s", title="pyr mem")
    p2=@df @subset(iso_sum_celltype_per, :interneuron .== true, :cuemem .== 1) scatter(:isolated_mean, :isolated_events_per_time, :events_per_time, group=:correct, alpha=0.5, camera=(i,30),ylabel="isolated multiunit/s", zlabel="multiunit/s", title="int mem")
    p3=@df @subset(iso_sum_celltype_per, :interneuron .== false, :cuemem .== 0) scatter(:isolated_mean, :isolated_events_per_time, :events_per_time, group=:correct, alpha=0.5, camera=(i,30), xlabel="fraction(isolated)",ylabel="isolated multiunit/s", zlabel="multiunit/s", title="pyr cue")
    p4=@df @subset(iso_sum_celltype_per, :interneuron .== true, :cuemem .== 0) scatter(:isolated_mean, :isolated_events_per_time, :events_per_time, group=:correct, alpha=0.5, camera=(i,30),ylabel="isolated multiunit/s", zlabel="multiunit/s", title="int cue")
    plot(p3,p4,p1,p2; layout=grid(2,2))
end
gif(anim, plotsdir(Plot.folder_args..., "separatrix-three-vars-gif_facet=pyrint,cuemem_group=correct.gif"); loop=1)

p1=@df @subset(iso_sum_celltype_per, :interneuron .== false, :cuemem .== 0) scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="pyr cue")
p2=@df @subset(iso_sum_celltype_per, :interneuron .== true, :cuemem .== 0) scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="int cue")
p3=@df @subset(iso_sum_celltype_per, :interneuron .== false, :cuemem .== 1) scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="pyr mem")
p4=@df @subset(iso_sum_celltype_per, :interneuron .== true, :cuemem .== 1) scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="pyr cue")
plot(p1,p2,p3,p4, markersize=2)
Plot.save((;desc="fract vs mua-per-sec"))
# ======================================================


"""
Mean cycle dist
"""
Plot.setfolder("nonlocality","isolation-basic_column_plots")
histogram(filter(x->!ismissing(x) && x<100, @subset(spikes,:area.=="CA1", :interneuron .!= true).meancyc), bins=40, yscale=:log10)
Plot.save((;desc="Mean cycle, all cells"))
histogram(filter(x->!ismissing(x) && x<100, @subset(spikes,:area.=="CA1", :interneuron .!= true).nearestcyc), bins=40, yscale=:log10)
Plot.save((;desc="Nearest cycle, all cells"))

"""
SVM: All three variables
"""
Plot.setfolder("nonlocality", "svm")
using MLJ
@load SVC pkg=LIBSVM
#cv = StratifiedCV(nfolds=2)
cv = Holdout(fraction_train=0.80, shuffle=true)

function eva(cv::ResamplingStrategy, svc::Machine; n=100)
    measure=[accuracy, measure_0, measure_1]
    if typeof(cv) <: StratifiedCV
        res = MLJ.evaluate!(svc; resampling=cv,
                           verbosity=0,measure)
        df = DataFrame(res)
    elseif typeof(cv) <: Holdout
        res = [MLJ.evaluate!(svc; resampling=cv,
                           verbosity=0, measure)
              for _ in 1:n]
        df = [transform!(DataFrame(r), :measure=>(x->i*ones(size(DataFrame(r),1)))=>:replicate) 
              for (i,r) in enumerate(res)]
        df = vcat(df...)
    end
    res, df
end
function measure_0(y,ŷ)
    inds = y .== 0
    accuracy(y[inds], ŷ[inds])
end
function measure_1(y,ŷ)
    inds = y .== 1
    accuracy(y[inds], ŷ[inds])
end
function get_confusion_matrix()
end
cmt=Dict("accuracy"=>"accuracy", "measure_0 "=>"pred. cue\npred. error", "measure_1 "=>"pred. memory\npred. correct")

svm_data = @subset(iso_sum_celltype_per, :correct .== 0 .|| :correct .== 1);
svm_data = dropmissing(svm_data[:,[:isolated_mean,:events_per_time,:isolated_events_per_time,:correct, :cuemem]]);
X = MLJ.table(Matrix(svm_data[!, [:isolated_mean,:events_per_time,:isolated_events_per_time]]));

y = categorical(svm_data[!,:correct]);
svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
svc = machine(svc_mdl, X, y);
res, df = eva(cv, svc)
ld = Utils.mlj.measureidxDict(df)
xtick = (collect(values(ld)), getindex.([cmt],collect(keys(ld))))
@df df scatter(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)), :measurement; ylim=(0,1), label="Predict correct",
               xtick)

y = categorical(svm_data[!,:cuemem]);
svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
svc = machine(svc_mdl, X, y);
res, df = eva(cv, svc)
xtick = (collect(values(ld)), getindex.([cmt],collect(keys(ld))))
@df df scatter!(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)), :measurement; ylim=(0,1), label="Predict cuemem",
               xtick)
Plot.save("all three vars")

"""
SVM: events_per_time
"""

svm_data = @subset(iso_sum_celltype_per, :correct .== 0 .|| :correct .== 1);
svm_data = dropmissing(svm_data[:,[:events_per_time,:correct, :cuemem]]);
X = MLJ.table(Matrix(svm_data[!, [:events_per_time,]]))

y = categorical(svm_data[!,:correct]);
svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
svc = machine(svc_mdl, X, y);
res, df = eva(cv, svc)
xtick = (collect(values(ld)), cmt.[collect(keys(ld))])
@df df scatter(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)), :measurement; ylim=(0,1), label="Predict correct",
               xtick)

y = categorical(svm_data[!,:cuemem]);
svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
svc = machine(svc_mdl, X, y);
res, df = eva(cv, svc)
xtick = (collect(values(ld)), getindex.([cmt],collect(keys(ld))))
@df df scatter!(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)), :measurement; ylim=(0,1), label="Predict cuemem",
               xtick)
Plot.save("events per time")

"""
SVM: isolated_mean
"""

svm_data = @subset(iso_sum_celltype_per, :correct .== 0 .|| :correct .== 1);
svm_data = dropmissing(svm_data[:,[:isolated_mean,:correct, :cuemem]]);
X = MLJ.table(Matrix(svm_data[!, [:isolated_mean,]]))

y = categorical(svm_data[!,:correct]);
svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
svc = machine(svc_mdl, X, y);
res, df = eva(cv, svc)
xtick = (collect(values(ld)), getindex.([cmt],collect(keys(ld))))
@df df scatter(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)), :measurement; ylim=(0,1), label="Predict correct",
               xtick)

y = categorical(svm_data[!,:cuemem]);
svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
svc = machine(svc_mdl, X, y);
res, df = eva(cv, svc)
xtick = (collect(values(ld)), getindex.([cmt],collect(keys(ld))))
@df df scatter!(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)), :measurement; ylim=(0,1), label="Predict cuemem",
               xtick)
Plot.save("isolation")


"""
Relation of isolated to out of field?

Answer: Confusing or not much! At least explored with hulls. Maybe still fits Jai's
"""
isoin_sum = combine(groupby(spikes, [:isolated,:interneuron]),:infield => x->mean(skipmissing(x)), renamecols=false)
dropmissing!(isoin_sum)

@df isoin_sum bar(:isolated, :infield, group=:interneuron; xticks=([0,1],["adj","iso"]), legend_title="IsInterneuron")

isoin_sum_celltype = combine(groupby(spikes, [:unit, :isolated, :interneuron]),:infield=>x->nanmean(skipmissing(x));renamecols=false)
isoin_sum_celltype[!,:infield] = replace(isoin_sum_celltype.infield, NaN=>missing)
dropmissing!(isoin_sum_celltype)

Plot.setfolder("nonlocality","isolation-infield")
@df @subset(isoin_sum_celltype,:interneuron .!= true) boxplot(:isolated, :infield, xticks=([0,1],["adj","iso"]))
@df @subset(isoin_sum_celltype,:interneuron .!= true) scatter!(:isolated + randn(size(:isolated)).*0.1, :infield, xticks=([0,1],["adj","iso"]))
Plot.save((;desc="Not much difference adjacent-iso in percent infield spikes",hullmax=1,thresh=0.8))


"""
Isolated spikes ISI
"""
Munge.spiking.nextandprev!(spikes)
perc = round(mean(spikes.neard .> .4)*100,sigdigits=2)


# =========PLOTS =======================================
Plot.setfolder("nonlocality", "isolated_ca1_rate_pfc")
histogram(log10.(spikes.neard),label="nearest spike")
vline!([log10(.400)], label="thresh", title="$perc % isolated")
Plot.save("nearest spike isolation")
# ======================================================


"""
PFC FIRING DURING CA1 ISOLATION VERSUS ADJACENT
"""
using Infiltrator

function get_pos_labels(group; grouping, labels)
    p = []
    grouping = !(grouping isa Vector) ? [grouping] : grouping
    for g in grouping
        push!(p, Symbol(g) => group[1, g])
    end
    for (k,lab) in labels
        push!(p, Symbol(String(k) * "_label") => lab[group[1, k]])
    end
    #:cuemem => group.cuemem[1],
    #:cuemem_label => clab[group.cuemem[1]],
    p
end

function get_pfcrate_samples_at_spike(spikes, R; grouping=[], labels=Dict())
    pfc_rate_isoAdjacent_meanOfTimes = DataFrame()
    spikes = grouping == [] ? spikes : dropmissing(spikes, grouping)
    @time @showprogress for group in groupby(spikes, [:isolated,grouping...])
        #@info "samples" size(group.time)
        sample = [R[time=B, unit=At(pfc_units)]
                    for B in Between.(group.time.-0.015, group.time.+0.015)]
        sample = [s for s in sample if !isempty(s)]
        sample = vcat(sample...)
        unit = repeat(vec(pfc_units)', size(sample,1))
        append!(pfc_rate_isoAdjacent_meanOfTimes,
        DataFrame(OrderedDict(
            :unit => vec(unit),
            :rate => vec(sample),
            get_pos_labels(group; grouping=[:isolated, grouping...], labels)...
           )))
    end
    #@time pfc_rate_isoAdjacent_meanOfTimes.rankrate = sortperm(pfc_rate_isoAdjacent_meanOfTimes.rate)
    pfc_rate_isoAdjacent_meanOfTimes
end

function compute_differences_df(pfc_rate_isoAdjacent; 
        grouping=[], 
        addsamps=false,
        labels=Dict())
    D = DataFrame()
    grouping = Symbol.(grouping)
    groups = groupby(pfc_rate_isoAdjacent, grouping)
    @info "groups" length(groups)
    for group in groups
         adjacent_spikes = @subset(group,:isolated.==0).rate
         iso_spikes      = @subset(group,:isolated.==1).rate
         #radjacent_spikes = @subset(group,:isolated.==0).rankrate
         #riso_spikes      = @subset(group,:isolated.==1).rankrate
         zadjacent_spikes = zscore(@subset(group,:isolated.==0).rate)
         ziso_spikes      = zscore(@subset(group,:isolated.==1).rate)
         test = missing
         pval = missing
         rtest = missing
         rpval = missing
         try
             test =  UnequalVarianceTTest(iso_spikes, adjacent_spikes)
             pval = pvalue(test,tail=:right)
             #rtest =  UnequalVarianceTTest(riso_spikes, radjacent_spikes)
             #rpval = pvalue(rtest,tail=:right)
         catch
            continue
         end
         sampgroups = addsamps ? [
            :samp_iso => [iso_spikes],
            :samp_adj => [adjacent_spikes],
           ] : []
         
         append!(D, OrderedDict(
            :diff            => mean(iso_spikes) - mean(adjacent_spikes),
            #:rankdiff        => median(riso_spikes) - median(radjacent_spikes),
            :zdiff           => mean(ziso_spikes) - mean(zadjacent_spikes),
            :zdiv            => mean(ziso_spikes)/mean(zadjacent_spikes),
            :div             => mean(iso_spikes)/mean(adjacent_spikes),
            :iso_spikes      => mean(iso_spikes),
            :adjacent_spikes => mean(adjacent_spikes),
            :pval => pval,
            :test => test,
            :rpval => rpval,
            :rtest => rtest,
            get_pos_labels(group; grouping, labels)...,
            sampgroups...
           ))
    end
    D
end

#function 

pfc_rate_isoAdjacent_meanOfTimes = get_pfcrate_samples_at_spike(spikes, R;
                                                                grouping=[],
                                                                labels=[])

#transform!(pfc_rate_isoAdjacent_meanOfTimes, DataFrames.All(), :isolated => (x->[isonames[xx] for xx in x]) => :isolated_label)
@subset!(pfc_rate_isoAdjacent_meanOfTimes, :isolated .== 0 .|| :isolated .== 1)
srt = UnequalVarianceTTest(@subset(pfc_rate_isoAdjacent_meanOfTimes,:isolated.==0).rate,
                           @subset(pfc_rate_isoAdjacent_meanOfTimes,:isolated.==1).rate)


# =========PLOTS =======================================
Plot.setfolder("nonlocality", "isolated_ca1_rate_pfc")
@df pfc_rate_isoAdjacent_meanOfTimes boxplot(:isolated, :rate, title="TTest=$(round(pvalue(srt),sigdigits=2)), with N=$(srt.df) PFC cells",
                         xticks=([0,1], collect(values(isonames))))
@df pfc_rate_isoAdjacent_meanOfTimes scatter!(:isolated, :rate, group=:isolated)
Plot.save((;save_kws...,desc="meanmean_pfc_cell_rate"))
# ======================================================

"""
Mean of rates, per cell x time
"""
D = compute_differences_df(pfc_rate_isoAdjacent_meanOfTimes; grouping=[], 
                           labels=Dict(:isolated=>isonames))

@subset!(D, :isolated .== 0 .|| :isolated .== 1)
transform!(D, DataFrames.All(), :isolated => (x->[isonames[xx] for xx in x]) => :isolated_label)
srt = UnequalVarianceTTest(@subset(D,:isolated.==0).rate,
                          @subset(D,:isolated.==1).rate)


# =========PLOTS =======================================
Plot.setfolder("nonlocality", "isolated_ca1_rate_pfc")
@df isoadj_cellandtime boxplot(:isolated, :rate, title="Firing rate samples per spike\nTTest=$(round(pvalue(srt),sigdigits=2)) difference of individual cells",
                         xticks=([0,1], collect(values(isonames))), outliers=false)
Plot.save((;save_kws...,desc="firing rate samples per spike"))
# ======================================================

"""
DO SAME split by cue and memory!

(STOPPED UPDATING HERE)
"""
#TODO

srt   = OrderedDict()
diffs = OrderedDict()
@showprogress for U in groupby(isoadj_cellandtime, :unit)
    push!(srt,
            U.unit[1] => UnequalVarianceTTest(@subset(U,:isolated.==0).rate,
                                              @subset(U,:isolated.==1).rate))
    push!(diffs,
          U.unit[1] => mean(@subset(U,:isolated.==0).rate) - mean(@subset(U,:isolated.==1).rate))

end

# =========PLOTS =======================================
Plot.setfolder("nonlocality", "isolated_ca1_rate_pfc")
pv_higheradjacent = pvalue.(values(srt)) .< 0.05/length(srt) .&& values(diffs) .> 0
pv_higherisolated  = pvalue.(values(srt)) .< 0.05/length(srt) .&& values(diffs) .< 0
pv_nonsig    = pvalue.(values(srt)) .> 0.05/length(srt)
res =[collect(values(diffs)) pvalue.(values(srt))]
bar([0, 1, 2], mean.([pv_nonsig, pv_higheradjacent, pv_higherisolated]),
    xticks=([0, 1, 2], ["not bonferroni\nsig", "cells sig adjacent", "cells sig isolated"]),
    ylabel="Percent", title="Nonlocality: PFC firing conditioned on CA1 field")

Plot.save((;save_kws...,desc="bonferroni sig per cell, firing rate samples per spike"))
# ======================================================

"""
# Split by cuemem
"""
Utils.filtreg.register(beh,spikes, on="time", transfer=["cuemem"])

pfc_rate_isoAdjacent_meanOfTimes = get_pfcrate_samples_at_spike(spikes, R; grouping=[:cuemem], labels=Dict(:cuemem=>clab))
D  = compute_differences_df(pfc_rate_isoAdjacent_meanOfTimes; grouping=:cuemem, labels=Dict(:cuemem=>clab))
D.clabel = getindex.([clab], D.cuemem)

transform!(D, Colon(), :pval => Utils.pfunc => :pval_label)
#xticks=([-1,0,1], collect(values(clab)
@df D bar(:clabel, :diff, group=:clabel)
markers = text.(D.pval_label, :white)
@df D annotate!(:clabel, :diff * 1.025, markers)
#transform!(Du, DataFrames.All(), :pval => Utils.pfunc => :pval_label)
Plot.save((;save_kws..., desc="diff of iso adjacent in each state"))

"""
# Split by cuemem and pfc unit
"""
pfc_rate_isoAdjacent_meanOfTimes = get_pfcrate_samples_at_spike(spikes, R; grouping=[:cuemem,:unit], labels=Dict(:cuemem=>clab))
Du = compute_differences_df(pfc_rate_isoAdjacent_meanOfTimes; grouping=[:unit,:cuemem], labels=Dict(:cuemem=>clab))
Du.sig = (Du.pval .< 0.05/21)
Du.rsig = (Du.rpval .< 0.05/21)
@df Du scatter(:cuemem + 0.1.*randn(size(:cuemem)), :diff, xticks=([-1,0,1], collect(values(clab))), group=:area, legend_position=:outerbottomright, ylabel="Δ(iso-adj) \$MUA_{pfc}\$ hz")
hline!([0])

# SHEER NUMBER OF ISOLATED SPIKES (ALSO computer atop this script)
if_counts = combine(groupby(Du, [:cuemem, :isolated]), nrow)
subset!(if_counts, :isolated => x->(!).(isnan.(x)))
if_perc = combine(groupby(if_counts, [:cuemem]), :isolated, :nrow => (x->x./(sum(x))) => :perc)
iflab = Dict(0 => "out", 1=> "in")
clab = Dict(-1 => "nontask", 0 => "cue", 1=> "mem")
transform!(if_perc, DataFrames.All(), [:cuemem, :isolated] => ((x,y)->(getindex.([clab], x) .* " " .* getindex.([iflab], y))) .=> :cuemem_inoutfield)
sort!(if_perc, :isolated)
@df subset(if_perc, :isolated => (x->x.!=false)) bar(:cuemem_inoutfield, :perc, group=:isolated)

Plot.save((;save_kws..., desc="cellcounts sig"))

"""
# PFC firign rates norm by 0-1
"""
Plot.setfolder("nonlocality", "isolated_ca1_pfc_rates")
Rpfc =  Munge.spiking.torate(@subset(allspikes, :area .== "PFC"), beh)
Dpfc = DataFrame([0;Rpfc], :auto)
Dpfc.cuemem = beh.cuemem
pfc_cell_means = Matrix(combine(groupby(Dpfc, "cuemem"), 1:21 .=> mean))[:,2:end]
pfc_cell_means./maximum(pfc_cell_means,dims=1)
heatmap(["Nontask","Cue","Memory"], 1:size(pfc_cell_means,2), pfc_cell_means./maximum(pfc_cell_means,dims=1))
Plot.save("Average PFC cell firing rates per task type")

#bar(["Nontask","Cue","Memory"], mean(pfc_cell_means./maximum(pfc_cell_means,dims=1), dims=2))

Plot.save("Average of PFC cell firing, mean(norm0)")
bar(["Nontask","Cue","Memory"], mean(pfc_cell_means, dims=2))
Plot.save("Average of PFC cell firing, mean(rate)")

"""
# PFC firing rates norm by zscore
"""
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

