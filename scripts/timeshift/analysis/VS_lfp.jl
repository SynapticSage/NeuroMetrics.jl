#=
===================================================
Script Name     : VS_lfp.jl
Author          : Ryan Young
Created         : 05-03-2023
Last Modified   : 05-03-2023
===================================================

# Description
This script plots local field potential properties against timeshifted
optimal 𝛕 of each cell.
=#

using GoalFetchAnalysis
using DataFrames, DataFramesMeta, ProgressMeter, StatsBase, Statistics, StatsPlots, Plots

#    _  _     _                    _       _       _        
#  _| || |_  | |    ___   __ _  __| |   __| | __ _| |_ __ _ 
# |_  ..  _| | |   / _ \ / _` |/ _` |  / _` |/ _` | __/ _` |
# |_      _| | |__| (_) | (_| | (_| | | (_| | (_| | || (_| |
#   |_||_|   |_____\___/ \__,_|\__,_|  \__,_|\__,_|\__\__,_|
#                                                           
animal, day = "super_clean", 0
CELLS  = DI.load_cells(animal, day, "*")
time_factors = DI.behavior_0time_before_process("super")
taus = names(CELLS)[occursin.("1of20", names(CELLS))]

# -------------------------------
# Copy pasta from setup_pyr_iso.jl
# -------------------------------
function get_data(animal, day)
    # l_pyr = DI.load_lfp(animal, day; append="pyr")
    # l_pyr = transform(l_pyr, :time=> t-> t .- time_factors[animal], renamecols=false)
    # transform!(l_pyr, :time=> t-> t .+ time_factors[animal], renamecols=false)
    cycles = DI.load_cycles(animal, day, "pyr")
    cycles.start = cycles.start_function
    cycles.stop  = cycles.stop_function
    cycles = transform(cycles,
        :start=> t-> t .- time_factors[animal],
        :stop=> t-> t .-  time_factors[animal],
        renamecols=false
    )
    cycles = cycles[!,Not([:start_function,:stop_function])]
    spikes = transform(DI.load_spikes(animal, day, "pyr_cycles_isolated"),
        :time=> t-> t .- time_factors[animal],
        renamecols=false
    )
    # spikes = subset(SPIKES, :animal=>a->a.==animal, :day=>d->d.==day)
    # Print out extrema(time) for each of these just to be sure
    begin
        println("spikes.time: ", extrema(spikes.time))
        # println("lfp.time: ", extrema(l_pyr.time))
        println("cycles.time: ", extrema(cycles.start))
    end
    DIutils.pushover("Loaded checkpoint for $animal")
    (;cycles, spikes)
end

CYCLES, SPIKES = [], [], []
datasets = [("RY16",36), ("RY22",21)]
(animal, day) = first(datasets)
for (animal, day) in datasets
    if !isfile(DI.cyclepath(animal, day, "pyr"))
        continue
    end
    println("Loading $animal")
    cycles, spikes = get_data(animal, day)
    println("Adding $taus to $animal")
    spikes = sort!(spikes, [:unit, :time])
    DIutils.filtreg.register(CELLS, spikes, on="unit", transfer=taus)
    println("Done with $animal")
    push!(CYCLES, cycles)
    push!(SPIKES, spikes)
    # Transfer the tau values to the spikes
end
CYCLES = vcat(CYCLES...)
SPIKES = vcat(SPIKES...)

#    _  _        _                _               
#  _| || |_     / \   _ __   __ _| |_   _ _______ 
# |_  ..  _|   / _ \ | '_ \ / _` | | | | |_  / _ \
# |_      _|  / ___ \| | | | (_| | | |_| |/ /  __/
#   |_||_|   /_/   \_\_| |_|\__,_|_|\__, /___\___|
#                                   |___/         

# First, let's actually create views into the spikes that actually occur
# in ripple and valid theta cycle
# ISSUE: better filter theta cycles before this step?
sptheta = subset(SPIKES, :theta=>t->t .!== missing .&& t.>0, view=true)
sprip   = subset(SPIKES, :ripple=>t->t .!== missing .&& t.>0,
        :ripple => r-> r.!= NaN, view=true)
spthrip = subset(SPIKES, :theta=>t->t .!== missing .&& t.>0, 
                         :ripple => r-> r.!= NaN, 
                         :ripple=>t->t .!== missing .&& t.>0, view=true)

gr()
histogram(filter(t->t!=0, sptheta.theta_phase))
histogram(filter(r->r!=0,sprip.ripple_phase))

histogram2d(spthrip.ripple_phase, spthrip.theta_phase, markersize=0.1, alpha=0.1, colorbar=true)

th = sort(dropmissing(combine(groupby(sptheta, :unit), :theta_phase=>mean,
    taus.=>first, renamecols=false), "all-1of20"),:theta_phase)
rip = sort(dropmissing(combine(groupby(sprip, :unit), :ripple_phase=>mean,
    taus.=>first, renamecols=false), "all-1of20"),:ripple_phase)
thrip = sort(dropmissing(combine(groupby(spthrip, :unit), :ripple_phase=>mean,
    :theta_phase=>mean, taus.=>first, renamecols=false), "all-1of20"),:ripple_phase)

scatter(th.theta_phase, th[!,"all-1of20"], markersize=0.1, alpha=0.1)
scatter(rip.ripple_phase, rip[!,"all-1of20"], markersize=0.1, alpha=0.1)
scatter(thrip.ripple_phase, thrip.theta_phase, thrip[!,"all-1of20"], markersize=0.1, alpha=0.1)
scatter(thrip.ripple_phase, thrip.theta_phase, markersize=0.1, alpha=0.1)


