
if !(:lfp in names(Main))
    include("../load_isolated.jl")
end

# cycles.time = (cycles.stop - cycles.start)/2 + cycles.start

lfp[!,:broadraw] = lfp[!,:smoothbroadraw]

#Check that cycle labels match in our 3 relvent data sources
begin
    using Blink, Interact
    colors = theme_palette(:auto)
    using Colors
    function getcol(x)
        if ismissing(x)
            RGBA(NaN, NaN, NaN, NaN)
        else
            colors[mod(x,length(colors))+1]
        end
    end

    w = Window()
    R = collect(zip(1:30:(size(cycles,1)-30), 30:30:(size(cycles,1))))
    ui = @manipulate for selector=eachindex(R)

        Nstart,Nstop = R[selector]
        cycs = cycles[Nstart:Nstop, :]
        spans = collect.(zip(cycs.start, cycs.stop))
        icycs = cycs.cycle
        p=plot(size=(2000,1000))
        sp = @subset(spikes, :cycle .>= minimum(icycs), :cycle .<= maximum(icycs))

        if !isempty(sp)
            @df sp begin
                scatter!(:time, :unit, c=getcol.(:cycle), markerstrokecolor=colorant"black")
            end
        end

        [(Plots.vspan!(span, legend=false, alpha=0.1, c=getcol.(icyc)); 
          Plots.annotate!(span[1], ylims()[mod(icyc,2)+1], text(icyc)))
         for (icyc, span) in zip(icycs, spans)]

        p

    end
    body!(w, ui)
end

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


#   _  _     ___           _       _           _ 
# _| || |_  |_ _|___  ___ | | __ _| |_ ___  __| |
#|_  ..  _|  | |/ __|/ _ \| |/ _` | __/ _ \/ _` |
#|_      _|  | |\__ \ (_) | | (_| | ||  __/ (_| |
#  |_||_|   |___|___/\___/|_|\__,_|\__\___|\__,_|
#                                                

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

# function plot_cycle()
# end

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
# pcs = [st.plotcycstat for st in STATS]
# res = [st.rangeerror for st in STATS]
# for (rs, pc) in zip(res, pcs)
#     display(plot(pc, background_color = rs ? :lightpink : :white))
# end



# --------------------------------------
#    _  _        _       _  _                      _    
#  _| || |_     / \   __| |(_) __ _  ___ ___ _ __ | |_  
# |_  ..  _|   / _ \ / _` || |/ _` |/ __/ _ \ '_ \| __| 
# |_      _|  / ___ \ (_| || | (_| | (_|  __/ | | | |_  
#   |_||_|   /_/   \_\__,_|/ |\__,_|\___\___|_| |_|\__| 
#                        |__/                           
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
