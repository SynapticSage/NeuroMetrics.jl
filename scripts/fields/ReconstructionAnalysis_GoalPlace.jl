
if dopoissonmodel
    compare = model.generate_comparison(F, ["goal","place"])
end

# PLACE-GOAL JOINT DISTRIBUTION P(X,Y,γ,p)
props = ["x", "y", "currentAngle", "currentPathLength"]
newkws = (; kws..., resolution=20, gaussian=0, props=props,
          filters=merge(kws.filters))
X = field.get_fields(beh, spikes; newkws...);
F["placegoal-joint"] = X

# Acquire marginals
F["place-marginal"] = operation.marginalize(X, dims=[3,4])
F["goal-marginal"]  = operation.marginalize(X, dims=[1,2])

