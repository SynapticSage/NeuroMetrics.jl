"""
Boilerplate for fields

plus runs all field analyses
"""
using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
# Grab our raw data
using DataFrames
using KernelDensity, Distributions
using Plots, Measures
using ProgressMeter
includet(srcdir("raw.jl"))
includet(srcdir("field.jl"))
includet(srcdir("filt.jl"))
spikes, beh, ripples, cells = raw.load("RY16", 36);
function sf(p, loc)
    p.attr[:size] = (1900, 1900)
    savefig(p, loc)
    savefig(p*".svg", loc)
end

# Generalized settings
resolution = 60; # field resolution
splitby=["unit", "area"]
kws=(;resolution, splitby, filters=merge(filt.speed, filt.cellcount))
runanalyses, ploton = false, false

if runanalyses
    include(scriptsdir("fields", "PlaceField.jl")) # 1
    include(scriptsdir("fields", "GoalField.jl")) # 2
    include(scriptsdir("fields", "GoalPlace.jl")) # 3
    include(scriptsdir("fields", "GoalPlace.jl")) # 4
end
