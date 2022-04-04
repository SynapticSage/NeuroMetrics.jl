includet(srcdir("operation.jl"))
includet(srcdir("model.jl"))

# PLACE-GOAL JOINT DISTRIBUTION P(X,Y,γ,p)
props = ["x", "y", "currentPathLength", "currentAngle"]
newkws = (; kws..., resolution=40, gaussian=0, props=props,
          filters=merge(kws.filters))
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

if ploton
    # Reconstructions
    field.plot.show_fields(R̂["place"])
    field.plot.show_fields(R̂["goal"])
    # Marginals
    field.plot.show_fields(F["place-marginal"].hist)
    field.plot.show_fields(F["goal-marginal"])
end
