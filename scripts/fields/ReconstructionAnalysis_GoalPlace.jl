includet(srcdir("model.jl"))
includet(srcdir("operation.jl"))

# PLACE-GOAL JOINT DISTRIBUTION P(X,Y,γ,p)
props = ["x", "y", "currentPathLength", "currentAngle"]
filters = merge(kws.filters,
                filt.correct,
                filt.notnan("currentAngle"), 
                filt.minmax("currentPathLength", 2, 150))
newkws = (; kws..., resolution=40, gaussian=0, props=props,
          filters=merge(kws.filters, filters))

X = field.get_fields(beh, spikes; newkws...);
F["placegoal-joint"] = X

# Acquire marginals P(X,Y), P(γ, p)
F["place-marginal"] = operation.marginalize(X, dims=[3,4])
F["goal-marginal"]  = operation.marginalize(X, dims=[1,2])
F["place-marginal-sq"] = operation.marginalize(X, dims=[3,4], dosqueeze=true)
F["goal-marginal-sq"]  = operation.marginalize(X, dims=[1,2], dosqueeze=true)

R̂ = Dict()
R̂["goal"] = operation.apply(model.reconstruction, F["placegoal-joint"].behdens, 
                            F["place-marginal"].hist)
R̂["place"] = operation.apply(model.reconstruction, F["placegoal-joint"].behdens, 
                            F["goal-marginal"].hist)

error = model.reconstruction_error(F["place-marginal-sq"], R̂["place"])

utils.pushover("Finished reconstruction")
if ploton
    # Reconstructions
    p = field.plot.show_fields(R̂["place"])
    savefig(p, plotsdir("fields", "reconstruction", "reconstructed_place_from_goal.svg"))
    savefig(p, plotsdir("fields", "reconstruction", "reconstructed_place_from_goal.png"))
    savefig(p, plotsdir("fields", "reconstruction", "reconstructed_place_from_goal.pdf"))
    p = field.plot.show_fields(R̂["goal"])
    savefig(p, plotsdir("fields", "reconstruction", "reconstructed_goal_from_place.svg"))
    savefig(p, plotsdir("fields", "reconstruction", "reconstructed_goal_from_place.png"))
    savefig(p, plotsdir("fields", "reconstruction", "reconstructed_goal_from_place.pdf"))
    # Marginals
    p = field.plot.show_fields(F["place-marginal-sq"].hist)
    savefig(p, plotsdir("fields", "reconstruction", "marginal_place.svg"))
    savefig(p, plotsdir("fields", "reconstruction", "marginal_place.png"))
    savefig(p, plotsdir("fields", "reconstruction", "marginal_place.pdf"))
    p = field.plot.show_fields(F["goal-marginal-sq"].hist)
    savefig(p, plotsdir("fields", "reconstruction", "marginal_goal.svg"))
    savefig(p, plotsdir("fields", "reconstruction", "marginal_goal.png"))
    savefig(p, plotsdir("fields", "reconstruction", "marginal_goal.pdf"))
end
