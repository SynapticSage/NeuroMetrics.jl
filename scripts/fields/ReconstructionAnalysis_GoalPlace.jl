quickactivate("/home/ryoung/Projects/goal-code/")
include(scriptsdir("fields", "Initialize.jl"))
using StatsPlots

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
X = operation.occnorm(X)

# Acquire marginals P(X,Y), P(γ, p)
F["place-marginal"] = operation.marginalize(X, dims=[3,4])
F["goal-marginal"]  = operation.marginalize(X, dims=[1,2])
F["place-marginal-sq"] = operation.marginalize(X, dims=[3,4], dosqueeze=true)
F["goal-marginal-sq"]  = operation.marginalize(X, dims=[1,2], dosqueeze=true)

# Obtain reconstructions!
R̂ = Dict()
R̂["goal"] = operation.apply(model.reconstruction, F["placegoal-joint"].occR, 
                            F["place-marginal"].Rₕ)
R̂["place"] = operation.apply(model.reconstruction, F["placegoal-joint"].occR, 
                            F["goal-marginal"].Rₕ)

# Get reconstruction model error summary
E_place_under_goal = model.reconstruction_error(F["place-marginal-sq"].Rₕ,
                                                R̂["place"])
E_goal_under_place = model.reconstruction_error(F["goal-marginal-sq"].Rₕ,
                                                R̂["goal"])
E_place_under_goal = table.to_dataframe(E_place_under_goal; name="error")
E_goal_under_place = table.to_dataframe(E_goal_under_place; name="error")
E = vcat(E_place_under_goal, E_goal_under_place; source=["pug","gup"])
E = vcat(E_place_under_goal, E_goal_under_place, source=:source=>["pug","gup"])
uE = unstack(E, :source, :error)
uE.ratio = uE.pug./uE.gup
uE.ratio_gt_1 = uE.ratio .>= 1

if ploton

    # Summaries
    @df E Plots.histogram(:error, group=:source)
    @df uE Plots.histogram(:ratio, label="Ratio of pug/gup", nbins=100)
    Plots.vline!([1], linestyle=:dash, c=:black, label="Model equivalence")
    @df uE Plots.histogram(:ratio_gt_1, label="Ratio of pug/gup", nbins=2, xticks=([0.25, 1.25],["Crappier under place","Crappier under goal"]))
    @df uE Plots.scatter(:pug, :gup, label="pug to gup", xlim=(0,3))
    Plots.plot!(0:30, 0:30, linestyle=:dash, c=:black, label="line of equivalence")

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
    p = field.plot.show_fields(F["place-marginal-sq"].Rₕ)
    savefig(p, plotsdir("fields", "reconstruction", "marginal_place.svg"))
    savefig(p, plotsdir("fields", "reconstruction", "marginal_place.png"))
    savefig(p, plotsdir("fields", "reconstruction", "marginal_place.pdf"))

    p = field.plot.show_fields(F["goal-marginal-sq"].Rₕ)
    savefig(p, plotsdir("fields", "reconstruction", "marginal_goal.svg"))
    savefig(p, plotsdir("fields", "reconstruction", "marginal_goal.png"))
    savefig(p, plotsdir("fields", "reconstruction", "marginal_goal.pdf"))


    # Missing samples
    Plots.heatmap(utils.squeeze(sum(Float64.(X.occzeroinds),dims=[3,4])), title="Quantification of missing samples\nin γ and p of (x,y,γ,p)\n")
    Plots.savefig(plotsdir("fields","reconstruction", "marginal_quantification_of_goal_missing_samples.svg"))

end
