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


clab = OrderedDict(-1 => "nontask", 0 => "cue", 1=> "mem", missing=>"sleep")
Munge.nonlocal.setclab(clab)
isonames =  OrderedDict(false => :adjacent, true=>:isolated)
filt_desc = OrderedDict(:all => "> 2cm/s")
save_kws = (;pfc_rate_analy=true)
filt = Filt.get_filters()
datacut = :all

animal, day = "RY16",36
Plot.setappend((;animal,day))
Plot.setparentfolder("nonlocality")

# ===================
# ACQUIRE DATA
# ===================
# Acquire data
@time spikes, beh, cells = Load.load(animal, day, data_source=["spikes","behavior", "cells"])
beh, spikes = Utils.filtreg.filterAndRegister(beh, spikes, on="time", transfer=["x","y","cuemem"], filters=filt[datacut], filter_skipmissingcols=true)
allspikes = copy(spikes)
beh2 = Load.load_behavior(animal,day)
Munge.nonlocal.setunfilteredbeh(beh2)

# Acquire LFP and isolated spikes
lfp = Load.load_lfp(animal, day, tet=5);
lfp.time = lfp.time .- Load.min_time_records[1]
lfp = Munge.lfp.annotate_cycles(lfp)
#sp = @subset(spikes, :tetrode .== 6);

F = load_fields()
kz = collect(filter(k->k.animal == animal && k.day == day, keys(F)))
@time f = F[bestpartialmatch(kz, 
                             (;datacut, widths=5,coactivity=nothing), 
                             nothing_means_removekey=true)];
f = ShiftedFields(deepcopy(f))
unitshift = Timeshift.types.matrixform(f)

# ===================
# OUT OF FIELD SPIKES
# ===================
# Setup a  shift-getting convenience method, the shifts, and a few metrics
shifts = collect(unitshift.dims[2])
push_dims!(unitshift)
push_celltable!( unitshift, cells, :unit, :area)
push_metric!( unitshift, Field.metrics.bitsperspike)
push_metric!( unitshift, Field.metrics.totalcount)
push_shiftmetric!(unitshift, best_tau!; metric=:bitsperspike)
annotate_nonlocal_spikes!(spikes, cells, unitshift, 0)

# ===================
# ISOLATED SPIKING
# ===================
#Munge.spiking.isolated(sp, lfp, include_samples=true)
Munge.spiking.isolated(spikes, lfp, include_samples=false)


# Which cells pass our criteria?
#region = :CA1
#metricfilter = metricfilters[region]
cycles = Munge.lfp.get_cycle_table(lfp)
Munge.lfp.annotate_cycles!(spikes, cycles)


"""
# ==========================================
# Visualize isolated spiking events per cell
# ==========================================
# """
# In this section, I'm visualizing events to make sure nothing is insance
#
# Ideal algo
#
# capture theta cycle stats
# - number of isolated spikes
# - distance to nearest cycles ahead and behind
#
# Visualize
# - theta
# - cycle cuts
# - spike raster

isospikes = @subset(spikes, :isolated, (!).(ismissing.(:isolated)),
                    (!).(ismissing.(:cycle)))
cycles.isounits   = Vector{Union{Missing,Vector}}(missing, size(cycles,1));
cycles.isotime    = Vector{Union{Missing,Vector}}(missing, size(cycles,1));
cycles.isoN       = Vector{Union{Missing,Int16}}(missing, size(cycles,1));
cycles.nearestcyc = Vector{Union{Missing,Vector}}(missing, size(cycles,1));
cycles.meancyc    = Vector{Union{Missing,Vector}}(missing, size(cycles,1));
cycles.cycLen     = Vector{Union{Missing, Int16}}(missing, size(cycles,1)); 
gs, gl, gc = groupby(isospikes,:cycle), groupby(lfp,:cycle), groupby(cycles,:cycle)
@showprogress for cycle in unique(disallowmissing(isospikes.cycle))
    s, l, c = gs[(;cycle)], gl[(;cycle)], gc[(;cycle)]
    c = view(c, 1, :)
    c[:isounits]   = s.unit
    c[:isotime]   = s.time
    c[:isoN]   = size(s,1)
    c[:nearestcyc] = disallowmissing(s.nearestcyc)
    c[:meancyc]    = disallowmissing(s.meancyc)
    c[:cycLen]     = size(l,1)
end


Plot.setfolder("examples, isolated cycles")
gs, gl, gc = groupby(isospikes,:cycle), groupby(lfp,:cycle), groupby(cycles,:cycle)
cycfilt = subset(dropmissing(cycles, :isounits),
                 :cycLen => l->l.>10)
cycle = first(eachrow(cycfilt))
cycleit = Iterators.take(Random.shuffle(eachrow(cycfilt)), 1_000)
@showprogress "plotting cycles" for cycle in cycleit
    stats = Plot.nonlocal.plot_cycle_example(gs,gl,cycle;return_stats=true)
    stats.plot
    Plot.save((;cycle=stats.cycle, units=stats.unit))
end



