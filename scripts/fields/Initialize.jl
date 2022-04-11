#
#Boilerplate for field
#
#plus runs all field analyses
#
using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
# Grab our raw data
using DataFrames
using KernelDensity, Distributions
using Plots, Measures
using ProgressMeter
using StatsPlots
using DataFramesMeta
using DataStructures: OrderedDict
includet(srcdir("raw.jl"))
includet(srcdir("field.jl"))
includet(srcdir("field/operation.jl"))
includet(srcdir("field/model.jl"))
includet(srcdir("filt.jl"))
includet(srcdir("table.jl"))
includet(srcdir("utils.jl"))
@time spikes, beh, ripples, cells = raw.load("RY16", 36);
function sf(p, loc)
    p.attr[:size] = (1900, 1900)
    savefig(p, loc)
    savefig(p*".svg", loc)
end

F = Dict() # Store field 
P = Dict() # Store poisson model
R̂ = Dict() # Store reconstructions

# Generalized settings
splitby=["unit", "area"]
kws=(;resolution, splitby, filters=merge(filt.speed_lib, filt.cellcount))
ploton, dofields, dopoissonmodel, doreconstruction = false, false, false, false

if dofields
    include(scriptsdir("fields", "PlaceField.jl")) # 1
    include(scriptsdir("fields", "GoalField.jl")) # 2
    include(scriptsdir("fields", "Poisson_CompareGoalPlace.jl")) # 3
    include(scriptsdir("fields", "PlotGoalPlace.jl")) # 4
    include(scriptsdir("fields", "SplitGoalPlace.jl")) # 5
end

if doreconstruction
    include(scriptsdir("fields", "ReconstructionAnalysis_GoalPlace_SplitBy.jl"))
end
