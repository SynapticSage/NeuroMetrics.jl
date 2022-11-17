quickactivate(expanduser("~/Projects/goal-code/")); using GoalFetchAnalysis
using Timeshift
using Timeshift.types
using Timeshift.shiftmetrics
using Field.metrics
using Plot
using Plot.receptivefield
using Utils.namedtup
using Munge.timeshift: getshift
using Munge.nonlocal
using Utils.statistic: pfunc
import Plot
using Munge.spiking
using Filt
using Infiltrator

using DataStructures: OrderedDict
using DimensionalData
import DimensionalData: Between
using ProgressMeter
using DataFrames, DataFramesMeta
using Statistics, NaNStatistics, StatsBase, StatsPlots, HypothesisTests, GLM
using Plots
using LazySets
using ElectronDisplay # to potentially defeat the GKS error

Utils.filtreg.register(beh,spikes;on="time",transfer=["velVec"])

"""
TODOS

- labels
    - titles
    - legend titles
    - any missing x/ylabels
- more accurate peak to peak
- split these by moving and still
- split by 1st and 2nd visit
"""


"""
===========================
Sheer diffferent rate of events (adjacent/isolated) over cue, mem, nontask
===========================
"""

iso_sum = get_isolation_summary(spikes)
sort!(iso_sum, [:area, :cuemem])

pfc_units = @subset(cells,:area.=="PFC").unit
R = Munge.spiking.torate(allspikes, beh)
pfc_units = intersect(pfc_units, collect(R.dims[2]))

isolated = last(groupby(subset(spikes, :isolated=>x->(!).(isnan.(x))) ,
                                        :isolated))
@assert all(isolated.isolated .== true)


Load.register(cells, spikes, on="unit", transfer=["meanrate"])
spikes.interneuron = spikes.meanrate .> 5
histogram(cells.meanrate)
iso_sum_celltype = get_isolation_summary(spikes,[:cuemem, :interneuron])
sort!(iso_sum_celltype, [:area, :interneuron, :cuemem])

Load.register(beh, spikes, on="time", transfer=["period","correct","velVec"])
@assert :period ∈ propertynames(spikes)
iso_sum_celltype_per = get_isolation_summary(spikes, [:cuemem, :interneuron,
                                                      :period, :correct])
sort!(iso_sum_celltype_per, [:area, :interneuron, :cuemem, :period])
@subset!(iso_sum_celltype_per, (:cuemem .== -1 .&& :correct .== -1) .||
                               (:cuemem .== 0 .&& :correct .!= -1)  .||
                               (:cuemem .== 1 .&& :correct .!= -1))
subset!(iso_sum_celltype_per, :events_per_time => x-> (!).(isinf.(x)))
iso_sum_celltype_per

# =========PLOTS =======================================
Plot.setfolder("MUA and isolated MUA")
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



# =========PLOTS =======================================
@df iso_sum_celltype bar(:cuearea, :events_per_time, ylabel="MUA events per second\n$(filt_desc[:all])", group=:interneuron; kws...)
Plot.save((;desc="MUA per second, celltype"))
@df iso_sum_celltype bar(:cuearea, :isolated_mean, ylabel="Isolated Spikes (sign of CA1-PFC interaction)\n$(filt_desc[:all])", group=:interneuron; kws...)
Plot.save((;desc="fraction of isolated spikes, celltype"))
@df iso_sum_celltype bar(:cuearea, :isolated_events_per_time, ylabel="Isolated MUA × sec⁻1\n(sign of CA1-PFC interaction)\n$(filt_desc[:all])", group=:interneuron; kws...)
Plot.save((;desc="isolated spikes per second, celltype"))
# ======================================================

"""
===========================
Period-wise calculations
===========================
"""

# =========PLOTS =======================================
# -------- mua per second ------------------------------
Plot.setfolder("MUA and isolated MUA", "period-wise")
@df @subset(iso_sum_celltype_per, :interneuron .== false) scatter(:cuemem .+ (randn(size(:cuemem)) .* 0.1), :events_per_time, group=:correct, xtick=([-1,0,1],["nontask","cue","mem"]),title="MUA events/time", xlabel="task", ylabel="events × time\$^{-1}\$, pyr cells", legend_title="correct/error", legend_position=:outerbottomright)
@df @subset(iso_sum_celltype_per, :interneuron .== true) scatter(:cuemem .+ (randn(size(:cuemem)) .* 0.1), :events_per_time, group=:correct, xtick=([-1,0,1],["nontask","cue","mem"]),title="MUA events/time", xlabel="task", ylabel="events × time\$^{-1}\$, pyr cells", legend_title="correct/error", legend_position=:outerbottomright)
Plot.save((;desc="mua-per-time",group=:correct,pyr=true))

@df @subset(iso_sum_celltype_per, :interneuron .== false) scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :events_per_time, group=:correct, alpha=0.5, xtick=([-1,0,1],["nontask","cue","mem"]),title="MUA events/time", xlabel="task", ylabel="events × time\$^{-1}\$, pyr cells", legend_title="correct/error", legend_position=:outerbottomright)
Plot.save((;desc="mua-per-time",group=:correct,pyr=true,sep_groups=true))

@df iso_sum_celltype_per scatter(:cuemem .+ (randn(size(:cuemem)) .* 0.1), :events_per_time, group=:interneuron, alpha=0.5, xtick=([-1,0,1],["nontask","cue","mem"]),title="MUA events/time", xlabel="task", ylabel="events × time\$^{-1}\$", legend_title="pyr/int", legend_position=:outerbottomright)
Plot.save((;desc="mua-per-time",group=:interneuron))

@df iso_sum_celltype_per scatter(:cuemem .+ 0.25.*(:interneuron .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :events_per_time, group=:interneuron, alpha=0.5, xtick=([-1,0,1],["nontask","cue","mem"]),title="MUA events/time", xlabel="task", ylabel="events × time\$^{-1}\$", legend_title="pyr/int", legend_position=:outerbottomright)
Plot.save((;desc="mua-per-time",group=:interneuron,group_sep=true))

@subset!(iso_sum_celltype_per, :events_per_time .!= Inf)

cpyr_ce = @df @subset(iso_sum_celltype_per, :area.=="CA1", :interneuron.==false) boxplot(:cuemem .+ 0.25.*(:correct .- 0.5), :events_per_time, group=:correct, alpha=0.5, title="ca1 pyr")
cint_ce = @df @subset(iso_sum_celltype_per, :area.=="CA1", :interneuron.==true) boxplot(:cuemem .+ 0.25.*(:correct .- 0.5), :events_per_time, group=:correct, alpha=0.5, title="ca1 int")

ppyr_ce = @df @subset(iso_sum_celltype_per, :area.=="PFC", :interneuron.==false) boxplot(:cuemem .+ 0.25.*(:correct .- 0.5), :events_per_time, group=:correct, alpha=0.5, title="pfc pyr")
#pint_ce = @df @subset(iso_sum_celltype_per, :area.=="PFC", :interneuron.==true) scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :isolated_mean, group=:correct, alpha=0.5, title="pfc int")
pint_ce = Plot.create_blank_plot()

plot(cpyr_ce, cint_ce, ppyr_ce, pint_ce, layout=grid(2,2))

# -------- percent isolated spikes ---------------------
@df @subset(iso_sum_celltype_per, :interneuron .== false) scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :isolated_mean, group=:correct, alpha=0.5, xtick=([-1,0,1],["nontask","cue","mem"]),title="Isolation", xlabel="task", ylabel="%percent isolated spikes", legend_title="correct/error", legend_position=:outerbottomright)

Plot.save((;desc="isolation", group=:correct, group_sep=true, ))

@df iso_sum_celltype_per scatter(:cuemem .+ 0.25.*(:interneuron .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :isolated_mean, group=:interneuron, alpha=0.5)

Plot.save((;desc="isolation", group=:interneuron, group_sep=true))

@df iso_sum_celltype_per scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :isolated_mean, group=:correct, alpha=0.5)

cpyr_ce = @df @subset(iso_sum_celltype_per, :area.=="CA1", :interneuron.==false) scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :isolated_mean, group=:correct, alpha=0.5, title="ca1 pyr")
cint_ce = @df @subset(iso_sum_celltype_per, :area.=="CA1", :interneuron.==true) scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :isolated_mean, group=:correct, alpha=0.5, title="ca1 int")

ppyr_ce = @df @subset(iso_sum_celltype_per, :area.=="PFC", :interneuron.==false) scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :isolated_mean, group=:correct, alpha=0.5, title="pfc pyr")
#pint_ce = @df @subset(iso_sum_celltype_per, :area.=="PFC", :interneuron.==true) scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :isolated_mean, group=:correct, alpha=0.5, title="pfc int")
pint_ce = Plot.create_blank_plot()

plot(cpyr_ce, cint_ce, ppyr_ce, pint_ce, layout=grid(2,2))
Plot.save((;desc="isolation cell type and correct"))

cpyr_ce = @df @subset(iso_sum_celltype_per, :area.=="CA1", :interneuron.==false) boxplot(:cuemem .+ 0.25.*(:correct .- 0.5), :isolated_mean, group=:correct, alpha=0.5, title="ca1 pyr")
cint_ce = @df @subset(iso_sum_celltype_per, :area.=="CA1", :interneuron.==true) boxplot(:cuemem .+ 0.25.*(:correct .- 0.5), :isolated_mean, group=:correct, alpha=0.5, title="ca1 int")

ppyr_ce = @df @subset(iso_sum_celltype_per, :area.=="PFC", :interneuron.==false) boxplot(:cuemem .+ 0.25.*(:correct .- 0.5), :isolated_mean, group=:correct, alpha=0.5, title="pfc pyr")
#pint_ce = @df @subset(iso_sum_celltype_per, :area.=="PFC", :interneuron.==true) scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :isolated_mean, group=:correct, alpha=0.5, title="pfc int")
pint_ce = Plot.create_blank_plot()

plot(cpyr_ce, cint_ce, ppyr_ce, pint_ce, layout=grid(2,2))
Plot.save((;desc="isolation cell type and correct"))


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
gif(anim, plotsdir(Plot.complete_folder_args..., "separatrix-three-vars-gif_pyr_facet=cuemem_group=correct.gif"); loop=1)

anim = @animate for i in 1:360
    p1=@df @subset(iso_sum_celltype_per, :interneuron .== false, :cuemem .== 1) scatter(:isolated_mean, :isolated_events_per_time, :events_per_time, group=:correct, alpha=0.5, camera=(i,30), xlabel="fraction(isolated)",ylabel="isolated multiunit/s", zlabel="multiunit/s", title="pyr mem")
    p2=@df @subset(iso_sum_celltype_per, :interneuron .== true, :cuemem .== 1) scatter(:isolated_mean, :isolated_events_per_time, :events_per_time, group=:correct, alpha=0.5, camera=(i,30),ylabel="isolated multiunit/s", zlabel="multiunit/s", title="int mem")
    p3=@df @subset(iso_sum_celltype_per, :interneuron .== false, :cuemem .== 0) scatter(:isolated_mean, :isolated_events_per_time, :events_per_time, group=:correct, alpha=0.5, camera=(i,30), xlabel="fraction(isolated)",ylabel="isolated multiunit/s", zlabel="multiunit/s", title="pyr cue")
    p4=@df @subset(iso_sum_celltype_per, :interneuron .== true, :cuemem .== 0) scatter(:isolated_mean, :isolated_events_per_time, :events_per_time, group=:correct, alpha=0.5, camera=(i,30),ylabel="isolated multiunit/s", zlabel="multiunit/s", title="int cue")
    plot(p3,p4,p1,p2; layout=grid(2,2))
end
gif(anim, plotsdir(Plot.complete_folder_args..., "separatrix-three-vars-gif_facet=pyrint,cuemem_group=correct.gif"); loop=1)

p1=@df @subset(iso_sum_celltype_per, :interneuron .== false, :cuemem .== 0, :area .== "PFC") scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="pyr cue")
p3=@df @subset(iso_sum_celltype_per, :interneuron .== false, :cuemem .== 1, :area .== "PFC") scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="pyr mem")
p2 = p4 = Plot.create_blank_plot()
plot(p1,p2,p3,p4, markersize=2)
Plot.save((;desc="fract vs mua-per-sec, PFC"))

p1=@df @subset(iso_sum_celltype_per, :interneuron .== false, :cuemem .== 0, :area .== "CA1") scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="pyr cue")
p2=@df @subset(iso_sum_celltype_per, :interneuron .== true, :cuemem .== 0, :area .== "CA1") scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="int cue")
p3=@df @subset(iso_sum_celltype_per, :interneuron .== false, :cuemem .== 1, :area .== "CA1") scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="pyr mem")
p4=@df @subset(iso_sum_celltype_per, :interneuron .== true, :cuemem .== 1, :area .== "CA1") scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="pyr cue")

plot(p1,p2,p3,p4, markersize=2)
Plot.save((;desc="fract vs mua-per-sec, CA1"))
# ======================================================

"""
===========================
summary of period iso_mean
===========================
"""
Plot.setfolder("MUA and isolated MUA", "period-wise-summary")


XX = combine(groupby(@subset(iso_sum_celltype_per, :interneuron .== false, :area .== "CA1"),
                     :cuemem), :isolated_mean => median, :isolated_mean => x->std(x)/sqrt(length(x))
            )
# TODO median, not mean, and bootstrap the median
kws=(;xlabel="cuemem", ylabel="isolated frac", title="CA1 pyr")
@df XX bar(:cuemem, :isolated_mean_median, yerror=:isolated_mean_function;
          kws...)
Plot.save(kws)

XX = combine(groupby(@subset(iso_sum_celltype_per, :interneuron .== true, :area .== "CA1"),
                     :cuemem), :isolated_mean => median, :isolated_mean => x->std(x)/sqrt(length(x))
            )
kws=(;xlabel="cuemem", ylabel="isolated frac", title="CA1 int")
@df XX bar(:cuemem, :isolated_mean_median, yerror=:isolated_mean_function;
          kws...)
Plot.save(kws)

XX = combine(groupby(@subset(iso_sum_celltype_per, :interneuron .== false, :area .== "PFC"),
                     :cuemem), :isolated_mean => median, :isolated_mean => x->std(x)/sqrt(length(x))
            )
kws=(;xlabel="cuemem", ylabel="isolated frac", title="PFC pyr")
@df XX bar(:cuemem, :isolated_mean_median, yerror=:isolated_mean_function;
          kws...)
Plot.save(kws)

"""
summary of period mua
"""
XXg = groupby(@subset(iso_sum_celltype_per, :area .== "CA1"), :cuemem)
XX = combine(XXg, :events_per_time => median, 
             :events_per_time => x->nanstd(x)./sqrt(length(x))
            )
# TODO median, not mean, and bootstrap the median
kws=(;xlabel="cuemem", ylabel="MUA per sec", title="CA1 MUA per second")
@df XX bar(:cuemem, :events_per_time_median, yerr=:events_per_time_function;
          linewidth=2, 
          kws...)
@df XX plot!(:cuemem, :events_per_time_median, yerr=:events_per_time_function, linestyle=Symbol("none"))
Plot.save(kws)

XXg = groupby(@subset(iso_sum_celltype_per, :interneuron .== false, :area .== "CA1"), :cuemem)
XX = combine(XXg, :events_per_time => median, 
             :events_per_time => x->nanstd(x)./sqrt(length(x))
            )
# TODO median, not mean, and bootstrap the median
kws=(;xlabel="cuemem", ylabel="MUA per sec", title="CA1 pyr MUA")
@df XX bar(:cuemem, :events_per_time_median, yerr=:events_per_time_function;
          linewidth=2, 
          kws...)
@df XX plot!(:cuemem, :events_per_time_median, yerr=:events_per_time_function, linestyle=:none)
Plot.save(kws)

XX = combine(groupby(@subset(iso_sum_celltype_per, :interneuron .== true, :area .== "CA1"),
                     :cuemem), 
             :events_per_time => median, 
             :events_per_time => x->std(x)/sqrt(length(x))
            )
kws=(;xlabel="cuemem", ylabel="MUA per sec", title="CA1 int MUA")
@df XX bar(:cuemem, :events_per_time_median, yerror=:events_per_time_function;
          kws...)
Plot.save(kws)

XX = combine(groupby(@subset(iso_sum_celltype_per, :interneuron .== false, :area .== "PFC"),
                     :cuemem), :events_per_time => median, :events_per_time => x->std(x)/sqrt(length(x))
            )
kws=(;xlabel="cuemem", ylabel="MUA per sec", title="PFC pyr")
@df XX bar(:cuemem, :events_per_time_median, yerror=:events_per_time_function;
          kws...)
Plot.save(kws)

"""
===========================
summary of iso events per time
===========================
"""

XX = combine(groupby(@subset(iso_sum_celltype_per, :interneuron .== false, :area .== "CA1"),
                     :cuemem), :isolated_events_per_time => median, :isolated_events_per_time => x->std(x)/sqrt(length(x))
            )
# TODO median, not mean, and bootstrap the median
kws=(;xlabel="cuemem", ylabel="iso events per sec", title="CA1 pyr")
@df XX bar(:cuemem, :isolated_events_per_time_median, yerror=:isolated_events_per_time_function;
          kws...)
Plot.save(kws)

XX = combine(groupby(@subset(iso_sum_celltype_per, :interneuron .== true, :area .== "CA1"),
                     :cuemem), :isolated_events_per_time => median, :isolated_events_per_time => x->std(x)/sqrt(length(x))
            )
kws=(;xlabel="cuemem", ylabel="iso events per sec", title="CA1 int")
@df XX bar(:cuemem, :isolated_events_per_time_median, yerror=:isolated_events_per_time_function;
          kws...)
Plot.save(kws)

XX = combine(groupby(@subset(iso_sum_celltype_per, :interneuron .== false, :area .== "PFC"),
                     :cuemem), :isolated_events_per_time => median, :isolated_events_per_time => x->std(x)/sqrt(length(x))
            )
kws=(;xlabel="cuemem", ylabel="iso events per sec", title="PFC pyr")
@df XX bar(:cuemem, :isolated_events_per_time_median, yerror=:isolated_events_per_time_function;
          kws...)
Plot.save(kws)


"""
Mean cycle dist
"""

Plot.setfolder("isolation-basic_column_plots")
histogram(filter(x->!ismissing(x) && x<100, @subset(spikes,:area.=="CA1", :interneuron .!= true).meancyc), bins=40, yscale=:log10)
Plot.save((;desc="Mean cycle, all cells"))
histogram(filter(x->!ismissing(x) && x<100, @subset(spikes,:area.=="CA1", :interneuron .!= true).nearestcyc), bins=40, yscale=:log10)
Plot.save((;desc="Nearest cycle, all cells"))

"""
===========================
SVM: All three variables
===========================
"""

Plot.setfolder( "svm")
using MLJ
MLJ.@load SVC pkg=LIBSVM
cv = StratifiedCV(nfolds=2, shuffle=true)
#cv = Holdout(fraction_train=0.80, shuffle=true)

function eva(cv::ResamplingStrategy, svc::Machine; n=100)
    measure=[accuracy, measure_0, measure_1]
    #if typeof(cv) <: StratifiedCV
    #    res = MLJ.evaluate!(svc; resampling=cv,
    #                       verbosity=0,measure)
    #    df = DataFrame(res)
    #elseif typeof(cv) <: Holdout
        res = [MLJ.evaluate!(svc; resampling=cv,
                           verbosity=0, measure)
              for _ in 1:n]
        df = [transform!(DataFrame(r), :measure=>(x->i*ones(size(DataFrame(r),1)))=>:replicate) 
              for (i,r) in enumerate(res)]
        df = vcat(df...)
    #end
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
XX = Matrix(svm_data[!, [:isolated_mean,:events_per_time,:isolated_events_per_time]]);

y = categorical(svm_data[!,:correct]);
svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
svc = machine(svc_mdl, XX, y);
res, df = eva(cv, svc)
ld = Utils.mlj.measureidxDict(df)
xtick = (collect(values(ld)), getindex.([cmt],collect(keys(ld))))
@df df scatter(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)), :measurement; ylim=(0,1), label="Predict correct",
               xtick)

y = categorical(svm_data[!,:cuemem]);
svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
svc = machine(svc_mdl, XX, y);
res, df = eva(cv, svc)
xtick = (collect(values(ld)), getindex.([cmt],collect(keys(ld))))
@df df scatter!(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)), :measurement; ylim=(0,1), label="Predict cuemem", legend=:none,
               xtick)


Plot.save("all three vars")

"""
===========================
SVM: events_per_time
===========================
"""

svm_data = @subset(iso_sum_celltype_per, :correct .== 0 .|| :correct .== 1);
svm_data = dropmissing(svm_data[:,[:events_per_time,:correct, :cuemem]]);
XX = MLJ.table(Matrix(svm_data[!, [:events_per_time,]]))

y = categorical(svm_data[!,:correct]);
svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
svc = machine(svc_mdl, XX, y);
res, df = eva(cv, svc)
xtick = (collect(values(ld)), getindex.([cmt], collect(keys(ld))))
@df df scatter(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)), :measurement; ylim=(0,1), label="Predict correct",
               xtick)

y = categorical(svm_data[!,:cuemem]);
svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
svc = machine(svc_mdl, XX, y);
res, df = eva(cv, svc)
@df df scatter!(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)), :measurement; ylim=(0,1), label="Predict cuemem",
               xtick)

Plot.save("events per time")

"""
===========================
SVM: isolated_mean
===========================
"""

svm_data = @subset(iso_sum_celltype_per, :correct .== 0 .|| :correct .== 1);
svm_data = dropmissing(svm_data[:,[:isolated_mean,:correct, :cuemem]]);
XX = MLJ.table(Matrix(svm_data[!, [:isolated_mean,]]))

y = categorical(svm_data[!,:correct]);
svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
svc = machine(svc_mdl, XX, y);
res, df = eva(cv, svc)
xtick = (collect(values(ld)), getindex.([cmt],collect(keys(ld))))
@df df scatter(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)), :measurement; ylim=(0,1), label="Predict correct",
               xtick)

y = categorical(svm_data[!,:cuemem]);
svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
svc = machine(svc_mdl, XX, y);
res, df = eva(cv, svc)
xtick = (collect(values(ld)), getindex.([cmt],collect(keys(ld))))
@df df scatter!(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)), :measurement; ylim=(0,1), label="Predict cuemem",
               xtick)

Plot.save("isolation")

"""
===========================
SVM: isolated_events_per_time
===========================
"""

svm_data = @subset(iso_sum_celltype_per, :correct .== 0 .|| :correct .== 1);
svm_data = dropmissing(svm_data[:,[:isolated_events_per_time,:correct, :cuemem]]);
XX = MLJ.table(Matrix(svm_data[!, [:isolated_events_per_time,]]))

y = categorical(svm_data[!,:correct]);
svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
svc = machine(svc_mdl, XX, y);
res, df = eva(cv, svc)
xtick = (collect(values(ld)), getindex.([cmt],collect(keys(ld))))
@df df scatter(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)), :measurement; ylim=(0,1), label="Predict correct",
               xtick)

y = categorical(svm_data[!,:cuemem]);
svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
svc = machine(svc_mdl, XX, y);
res, df = eva(cv, svc)
xtick = (collect(values(ld)), getindex.([cmt],collect(keys(ld))))
@df df scatter!(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)), :measurement; ylim=(0,1), label="Predict cuemem",
               xtick)
Plot.save("isolated_events_per_time")


"""
===========================
Relation of isolated to out of field?
===========================

Answer: Confusing or not much! At least explored with hulls. Maybe still fits Jai's
"""

isoin_sum = combine(groupby(spikes, [:isolated,:interneuron]),:infield => x->mean(skipmissing(x)), renamecols=false)
dropmissing!(isoin_sum)

@df isoin_sum bar(:isolated, :infield, group=:interneuron; xticks=([0,1],["adj","iso"]), legend_title="IsInterneuron")

isoin_sum_celltype = combine(groupby(spikes, [:unit, :isolated, :interneuron]),:infield=>x->nanmean(skipmissing(x));renamecols=false)
isoin_sum_celltype[!,:infield] = replace(isoin_sum_celltype.infield, NaN=>missing)
dropmissing!(isoin_sum_celltype)

Plot.setfolder("isolation-infield")
@df @subset(isoin_sum_celltype,:interneuron .!= true) boxplot(:isolated, :infield, xticks=([0,1],["adj","iso"]))
@df @subset(isoin_sum_celltype,:interneuron .!= true) scatter!(:isolated + randn(size(:isolated)).*0.1, :infield, xticks=([0,1],["adj","iso"]))
Plot.save((;desc="Not much difference adjacent-iso in percent infield spikes",hullmax=1,thresh=0.8))


"""
===========================
Isolated spikes ISI
===========================
"""

Munge.spiking.nextandprev!(spikes)
perc = round(mean(spikes.neard .> .4)*100,sigdigits=2)


# =========PLOTS =======================================
Plot.setfolder( "isolated_ca1_rate_pfc")
histogram(log10.(spikes.neard),label="nearest spike")
vline!([log10(.400)], label="thresh", title="$perc % isolated")
Plot.save("nearest spike isolation")
# ======================================================


