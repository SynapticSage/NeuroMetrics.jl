"""
Boilerplate for field

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
includet(srcdir("table.jl"))
@time spikes, beh, ripples, cells = raw.load("RY16", 36);
function sf(p, loc)
    p.attr[:size] = (1900, 1900)
    savefig(p, loc)
    savefig(p*".svg", loc)
end
F = Dict() # Store field 
P = Dict() # Store poisson model

# Generalized settings
resolution = 80; # field resolution
splitby=["unit", "area"]
kws=(;resolution, splitby, filters=merge(filt.speed_lib, filt.cellcount))
runanalyses, ploton, dopoissonmodel = true, false, true

if runanalyses
    include(scriptsdir("fields", "PlaceField.jl")) # 1
    include(scriptsdir("fields", "GoalField.jl")) # 2
    include(scriptsdir("fields", "Poisson_CompareGoalPlace.jl")) # 3
    include(scriptsdir("fields", "PlotGoalPlace.jl")) # 4
    include(scriptsdir("fields", "SplitGoalPlace.jl")) # 5
end
