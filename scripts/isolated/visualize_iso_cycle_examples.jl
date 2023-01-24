quickactivate(expanduser("~/Projects/goal-code/")); using GoalFetchAnalysis
using Infiltrator
using DataStructures: OrderedDict
using DimensionalData
import DimensionalData: Between
using ProgressMeter
using DataFrames, DataFramesMeta
using Statistics, NaNStatistics, StatsBase, StatsPlots, HypothesisTests, GLM
using Plots
using LazySets

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

Plot.off()


cycles.time = (cycles.stop - cycles.start)/2 + cycles.start
Utils.filtreg.register(cycles, spikes; on="time", transfer=["cycle"])


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

isospikes = @subset(spikes, 
                    :isolated, 
                    (!).(ismissing.(:isolated)),
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
    c[:isotime]    = s.time
    c[:isoN]       = size(s,1)
    c[:nearestcyc] = disallowmissing(s.nearestcyc)
    c[:meancyc]    = disallowmissing(s.meancyc)
    c[:cycLen]     = size(l,1)
end
dropmissing(cycles, :isounits)

function plot_cycle()
end

# --------------------------------
# ACTUAL ISOLATED SPIKING EXAMPLES
# --------------------------------
Plot.setfolder("examples, isolated cycles")
begin
    gs, gl, gc, gS = groupby(isospikes,:cycle), groupby(lfp,:cycle), groupby(cycles,:cycle), groupby(spikes,:cycle)

    cycfilt = subset(dropmissing(cycles, :isounits),
                     :cycLen => l->l.>10)
    @assert !isempty(cycfilt)

    cycle = first(eachrow(cycfilt))
    cycleit = Iterators.take(Random.shuffle(eachrow(cycfilt)), 100)
    STATS, pushstats = [], true
    @showprogress "plotting cycles" for cycle in cycleit
        stats = Plot.nonlocal.plot_cycle_example(gs,gl,cycle;
                                                 return_stats=true, gS)
        pushstats ? push!(STATS, stats) : nothing
        Plot.save((;cycle=stats.cycle, units=stats.unit))
    end
end
pcs = [st.plotcycstat for st in STATS]
res = [st.rangeerror for st in STATS]
for (rs, pc) in zip(res, pcs)
    display(plot(pc, background_color = rs ? :lightpink : :white))
end
# --------------------------------------

# --------------------------------------
Plot.setfolder("examples, adjacent cycles, indiv")
begin
    adjspikes = @subset(spikes, (!).(:isolated), (!).(ismissing.(:isolated)),
                        (!).(ismissing.(:cycle)))
    cycles.adjunits   = Vector{Union{Missing,Vector}}(missing, size(cycles,1));
    cycles.adjtime    = Vector{Union{Missing,Vector}}(missing, size(cycles,1));
    cycles.adjN       = Vector{Union{Missing,Int16}}(missing,  size(cycles,1));
    gs, gl, gc = groupby(adjspikes,:cycle), groupby(lfp,:cycle), groupby(cycles,:cycle)
    @showprogress for cycle in unique(disallowmissing(adjspikes.cycle))
        s, l, c = gs[(;cycle)], gl[(;cycle)], gc[(;cycle)]
        c = view(c, 1, :)
        c[:adjunits]   = s.unit
        c[:adjtime]   = s.time
        c[:adjN]   = size(s,1)
        c[:nearestcyc] = disallowmissing(s.nearestcyc)
        c[:meancyc]    = disallowmissing(s.meancyc)
        c[:cycLen]     = size(l,1)
    end

    gs, gl, gc = groupby(adjspikes,:cycle), groupby(lfp,:cycle), groupby(cycles,:cycle)
    cycfilt = subset(dropmissing(cycles, :adjunits))
    cycle   = first(eachrow(cycfilt))
    cycleit = Iterators.take(Random.shuffle(eachrow(cycfilt)), 100)
    @showprogress "plotting cycles" for cycle in cycleit
        stats = Plot.nonlocal.plot_cycle_example(gs,gl,cycle;return_stats=true, gS, )
        Plot.save((;cycle=stats.cycle, units=stats.unit))
    end
end
