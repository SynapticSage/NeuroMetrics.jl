quickactivate("/home/ryoung/Projects/goal-code/")
include(scriptsdir("fields", "Initialize.jl"))

using StatsPlots

# PLACE-GOAL JOINT DISTRIBUTION P(X,Y,γ,p)
props = ["x", "y", "currentPathLength", "currentAngle","stopWell"]
filters = merge(kws.filters,
                filt.correct,
                filt.notnan("currentAngle"), 
                filt.minmax("currentPathLength", 2, 150))
newkws = (; kws..., resolution=[40, 40, 40, 40], gaussian=0, props=props,
          filters=merge(kws.filters, filters))
X = field.get_fields(beh, spikes; newkws...);
F["placegoal-joint"] = X
X = operation.occnorm(X)
goal_dims = [3,4]
place_dims = [1,2]

# Acquire marginals P(X,Y), P(γ, p)
F["place-marginal"]    = operation.marginalize(X, dims=goal_dims)
F["goal-marginal"]     = operation.marginalize(X, dims=place_dims)
F["place-marginal-sq"] = operation.marginalize(X, dims=goal_dims, dosqueeze=true)
F["goal-marginal-sq"]  = operation.marginalize(X, dims=place_dims, dosqueeze=true)

# Obtain reconstructions!
R̂ = Dict()
R̂["goal"] = operation.apply(model.reconstruction,
                                  F["placegoal-joint"].occR,
                                  F["place-marginal"].Rₕ)
R̂["place"] = operation.apply(model.reconstruction,
                            F["placegoal-joint"].occR, 
                            F["goal-marginal"].Rₕ)

# Get reconstruction model error summary
E_place_under_goal = model.reconstruction_error(F["place-marginal-sq"].Rₕ, R̂["place"])
E_goal_under_place = model.reconstruction_error(F["goal-marginal-sq"].Rₕ, R̂["goal"])
E_place_under_goal = table.to_dataframe(E_place_under_goal; name="error"x
E_goal_under_place = table.to_dataframe(E_goal_under_place; name="error")
E = vcat(E_place_under_goal, E_goal_under_place; source=["pug","gup"])
E = vcat(E_place_under_goal, E_goal_under_place, source=:source=>["pug","gup"])
uE = unstack(E, :source, :error)
uE.PG_GP_ratio = uE.pug./uE.gup
uE.PG_GP_diff = uE.pug.-uE.gup
uE.PG_GP_ratio_gt1 = uE.ratio .>= 1
uE.PG_GP_ratio_gt1 = replace(uE.ratio .>= 1)
cells = leftjoin(cells, uE[:,[:unit,:pug,:gup,:PG_GP_ratio, :PG_GP_ratio_gt1]])


if ploton


    p1=@df @subset(E,:area.=="CA1") Plots.histogram(:error, group=:source,
                                                    alpha=.6, title="CA1",
                                                    nbins=100)
    p2=@df @subset(E,:area.=="PFC") Plots.histogram(:error, group=:source,
                                                    alpha=.6, title="PFC",
                                                    nbins=100, xlim=(0,0.65))
    Plots.plot(p1,p2)

    p1=@df @subset(uE,:area.=="CA1") Plots.histogram(:PG_GP_diff,
                                                     label="Difference of
                                                     pug/gup", nbins=100,
                                                     c=:lightblue, title="CA1",
                                                     xlabel="ε(p|g) - ε(g|p)")
    vline!([0],linestyle=:dash,c=:violet,linewidth=2, label="equivalence")
    p2=@df @subset(uE,:area.=="PFC") Plots.histogram(:PG_GP_diff,
                                                     label="Difference of
                                                     pug/gup", nbins=100,
                                                     c=:red, title="PFC",
                                                     xlabel="ε(p|g) - ε(g|p)")
    vline!([0],linestyle=:dash,c=:violet,linewidth=2, label="equivalence")
    Plots.plot(p1,p2)

    @df @subset(uE, :area.=="CA1") Plots.histogram(:PG_GP_ratio_gt1, label="Ratio of pug/gup", nbins=2)
    @df @subset(uE, :area.=="CA1") Plots.histogram(:PG_GP_ratio_gt1, label="Ratio of pug/gup", nbins=2)
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
    Plots.heatmap(utils.squeeze(sum(Float64.(X.occzeroinds),dims=goal_dims)), title="Quantification of missing samples\nin γ and p of (x,y,γ,p)\n")
    Plots.savefig(plotsdir("fields","reconstruction", "goal_marginalize_quantification_of_goal_missing_samples.svg"))
    Plots.heatmap(utils.squeeze(mean(Float64.(X.occzeroinds),dims=goal_dims)), title="Quantification of missing samples\nin γ and p of (x,y,γ,p)\n")
    Plots.savefig(plotsdir("fields","reconstruction", "goal_marginalize_FRACTION_quantification_of_goal_missing_samples.svg"))
    Plots.heatmap(utils.squeeze(sum(Float64.(X.occzeroinds),dims=place_dims)), title="Quantification of missing samples\nin x and y of (x,y,γ,p)\n")
    Plots.savefig(plotsdir("fields","reconstruction", "place_marginalize_quantification_of_goal_missing_samples.svg"))
    Plots.heatmap(utils.squeeze(mean(Float64.(X.occzeroinds),dims=place_dims)), title="Quantification of missing samples\nin x and y of (x,y,γ,p)\n")
    Plots.savefig(plotsdir("fields","reconstruction", "place_marginalize_FRACTION_quantification_of_goal_missing_samples.svg"))

end
