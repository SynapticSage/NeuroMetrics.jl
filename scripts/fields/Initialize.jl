#
#Boilerplate for field
#
#plus runs all field analyses
#
using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
include(scriptsdir("fields","Include.jl"))
@assert Field.get_fields isa Function
@time spikes, beh, ripples, cells = Load.load("RY16", 36);
function sf(p, loc)
    p.attr[:size] = (1900, 1900)
    savefig(p, loc)
    savefig(p*".svg", loc)
end

F = Dict() # Store field 
P = Dict() # Store poisson model
R̂ = Dict() # Store reconstructions

# Generalized settings
Filt    = GoalFetchAnalysis.Filt
splitby = ["unit", "area"]
kws=(;resolution=80, splitby, filters=Filt.get_filters()[:all])
ploton, dofields, dopoissonmodel, doreconstruction = false, false, false, false

if dofields
    include(scriptsdir("fields", "PlaceField.jl"))               # 1
    include(scriptsdir("fields", "GoalField.jl"))                # 2
    include(scriptsdir("fields", "Poisson_CompareGoalPlace.jl")) # 3
    include(scriptsdir("fields", "PlotGoalPlace.jl"))            # 4
    include(scriptsdir("fields", "SplitGoalPlace.jl"))           # 5
end

if doreconstruction
    include(scriptsdir("fields", "ReconstructionAnalysis_GoalPlace_SplitBy.jl"))
end
