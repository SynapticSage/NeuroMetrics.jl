begin
    quickactivate(expanduser("~/Projects/goal-code/")); 
    using GoalFetchAnalysis
    using Infiltrator, DimensionalData, ProgressMeter, DataFrames, DataFramesMeta,
          Statistics, NaNStatistics, StatsBase, StatsPlots, HypothesisTests, GLM, Plots,
          LazySets, JLD2
    import DimensionalData: Between
    using DataStructures: OrderedDict

    using Timeshift, Timeshift.types, Timeshift.shiftmetrics, Field.metrics, Plot,
          Plot.receptivefield, DIutils.namedtup, Munge.nonlocal, Munge.spiking, Filt
          
    import GoalFetchAnalysis.Munge.isolated: parser
    import Plot
    using Munge.timeshift: getshift
    using DIutils.statistic: pfunc
    Plot.off()
    opt = parser()
end



jldopen(path_iso(opt),"r") do storage
    DIutils.dict.load_dict_to_module!(Main, Dict(k=>storage[k] for k in keys(storage)))
end

cycles.time = (cycles.stop - cycles.start)/2 + cycles.start
DIutils.filtreg.register(cycles, spikes; on="time", transfer=["cycle"])

"""
Check that cycle labels match in our 3 relvent data sources
"""
using Blink, Interact
colors = theme_palette(:auto)

    using Colors
@manipulate for Nstart=1:30:(size(cycles,1)-30), Nstop=30:30:(size(cycles,1))
    cycs = (collect.(zip(cycles.start, cycles.stop))[1:Nstart])
    icycs = cycles.cycle[1:Nstart]
    p=plot()
    Plots.vspan!(cycs, legend=false, alpha=0.1)
    @df @subset(spikes, :cycle .>= 1, :cycle .<= Nstart) begin
        scatter!(:time, :unit, c=colors[:cycle], markerstrokecolor=colorant"black")
    end
    p
end


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
