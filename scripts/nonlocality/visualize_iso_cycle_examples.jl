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

