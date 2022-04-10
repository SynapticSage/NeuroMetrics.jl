quickactivate("/home/ryoung/Projects/goal-code/")
@time include(scriptsdir("fields", "Initialize.jl"))

# ----------------------------------------
# FUNCTION SPECIFIC SHORTCUTS AND SETTINGS
# ----------------------------------------
# Marginals and reconstructions, and data structures to map em out
shortcut_names = OrderedDict("x-y" => "place",
                     "currentAngle"=>"γ",
                     "currentPathLength"=>"p",
                     "stopWell"=>"G")
𝕀(d) = Dict(zip(values(shortcut_names), keys(shortcut_names)))[d]
reconstructions = (("place","goal"),
                   ("goal", "place"),
                   ("place","particular-goal"),
                   ("particular-goal", "place"),
                   ("γ", "place"),
                   ("aparticular-γ", "place"))
props = ["x", "y", "currentPathLength", "currentAngle","stopWell"]

# Shortcut functions
"""
Returns the integer dim indices for a prop-string e.g. "x-y"->[1,2]
"""
𝔻(dims) = [findfirst(dim.==props) for dim in split(dims,"-")]
"""
Returns the remaining dimensions not covered by a prop-string
"""
𝔻₀(dims) = setdiff(1:length(props), 𝔻(dims)) # dims inverse
"""
Returns remaning dimensions as a joined prop string, instead of ints
"""
𝔻₀ⱼ(dims) = join(props[𝔻₀(dims)],"-") # joined
"""
Returns string with shortcut names
and the ₀ version returns the remaining names
"""
ℝ(dims) = replace(dims, [shortcut_names...][begin:1:end]...) # replace
ℝ₀ⱼ(dims) = ℝ(join(props[𝔻₀(dims)],"-")) # joined


# ----------------------------------------
# PLACE-GOAL JOINT DISTRIBUTION P(x,y, γ,p,G)
# ----------------------------------------
filters = merge(kws.filters,
                filt.correct,
                filt.notnan("currentAngle"), 
                filt.minmax("currentPathLength", 2, 150))
newkws = (; kws..., resolution=[40, 40, 40, 40, 5], gaussian=0, props=props,
          filters=merge(kws.filters, filters))
@time X = field.get_fields(beh, spikes; newkws...);
F["placegoal-joint"] = X
X = operation.occnorm(X)

# ---------
# MARGINALS
# ---------
# Acquire marginals P(X,Y), P(γ, p, G)
@time @showprogress for (dims, name) in shortcut_names
    @time F["$name-marginal"]    = operation.marginalize(X, dims=𝔻₀(dims));
end

# ---------------
# Reconstructions
# ---------------
# Obtain reconstructions!
R̂ = Dict()
for (reconstruction, marginal) in reconstructions
    marginal = ℝ(𝔻₀ⱼ(𝕀(reconstruction)))
    @time R̂[reconstruction] = operation.apply(model.reconstruction,
                                          F["placegoal-joint"].occR,
                                          F["$marginal-marginal"].Rₕ);
end

# ---------------
# SUMMARIES
# ---------------
# Get reconstruction model error summary
E_place_under_goal = model.reconstruction_error(F["place-marginal"].Rₕsq, R̂["place"])
E_goal_under_place = model.reconstruction_error(F["goal-marginal"].Rₕsq,  R̂["goal"])
E_place_under_goal = table.to_dataframe(E_place_under_goal; name="error")
E_goal_under_place = table.to_dataframe(E_goal_under_place; name="error")
E = vcat(E_place_under_goal, E_goal_under_place; source=["pug","gup"])
E = vcat(E_place_under_goal, E_goal_under_place, source=:source=>["pug","gup"])
uE = unstack(E, :source, :error)
uE.PG_GP_ratio = uE.pug./uE.gup
uE.PG_GP_diff = uE.pug.-uE.gup
uE.PG_GP_ratio_gt1 = uE.PG_GP_ratio .>= 1
uE.PG_GP_ratio_gt1_str = replace(uE.PG_GP_ratio .>= 1)
cells = leftjoin(cells, uE[:,[:unit,:pug,:gup,:PG_GP_ratio, :PG_GP_ratio_gt1]])


                    
#,---.|         |    
#|---'|    ,---.|--- 
#|    |    |   ||    
#`    `---'`---'`---'
if ploton


    ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
    #           SUMMARIES ... of reconstructions ...       #
    ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
    
    # Title: Visaluze 𝓍 = error 𝒷𝓎 {PUG, GUP} x (AREA)
    p1=@df @subset(E,:area.=="CA1") Plots.histogram(:error, group=:source,
                                                    alpha=.6, title="CA1",
                                                    nbins=100)
    p2=@df @subset(E,:area.=="PFC") Plots.histogram(:error, group=:source,
                                                    alpha=.6, title="PFC",
                                                    nbins=100, xlim=(0,0.65))
    Plots.plot(p1,p2)

    # Title: Visaluze 𝓍 = Δ(pug-gup) 𝒷𝓎 (AREA)
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


    ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
    #              RECONSTRUCTIONS                         #
    ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
    
    p = field.plot.show_fields(R̂["place"])
    savefig(p, plotsdir("fields", "reconstruction", "reconstructed_place_from_goal.svg"))
    savefig(p, plotsdir("fields", "reconstruction", "reconstructed_place_from_goal.png"))
    savefig(p, plotsdir("fields", "reconstruction", "reconstructed_place_from_goal.pdf"))
    p = field.plot.show_fields(R̂["goal"])
    savefig(p, plotsdir("fields", "reconstruction", "reconstructed_goal_from_place.svg"))
    savefig(p, plotsdir("fields", "reconstruction", "reconstructed_goal_from_place.png"))
    savefig(p, plotsdir("fields", "reconstruction", "reconstructed_goal_from_place.pdf"))


    ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
    #                   MARGNIAL                           #
    ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##

    p = field.plot.show_fields(F["place-marginal"].Rₕsq)
    savefig(p, plotsdir("fields", "reconstruction", "marginal_place.svg"))
    savefig(p, plotsdir("fields", "reconstruction", "marginal_place.png"))
    savefig(p, plotsdir("fields", "reconstruction", "marginal_place.pdf"))
    p = field.plot.show_fields(F["goal-marginal"].Rₕsq)
    savefig(p, plotsdir("fields", "reconstruction", "marginal_goal.svg"))
    savefig(p, plotsdir("fields", "reconstruction", "marginal_goal.png"))
    savefig(p, plotsdir("fields", "reconstruction", "marginal_goal.pdf"))

    ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
    #     JUXTAPOSITION OF RECONSTRUCTION AND MARGNALS     #
    ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##

    groups=field.group([F["place-marginal"].Rₕsq, R̂["place"]], 
                       ["M(place)", "R̂(place)"])
    overall = field.plot.show_fieldgroups(groups)

    # SEGFAULTS w/o filtering
    #groups=field.group([F["goal-marginal-sq"].Rₕ, R̂["goal"]], 
    #                   ["M(place)", "R̂(place)"])
    #overall = field.plot.show_fieldgroups(groups)
    


    ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
    #                   SUMMARY OF MISSING SAMPLES         #
    ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
    goal_dims = [3,4,5]
    place_dims = [1,2]
    Plots.heatmap(utils.squeeze(sum(Float64.(X.occzeroinds),dims=goal_dims)), title="Quantification of missing samples\nin γ and p of (x,y,γ,p)\n")
    Plots.savefig(plotsdir("fields","reconstruction", "goal_marginalize_quantification_of_goal_missing_samples.svg"))
    Plots.heatmap(utils.squeeze(mean(Float64.(X.occzeroinds),dims=goal_dims)), title="Quantification of missing samples\nin γ and p of (x,y,γ,p)\n")
    Plots.savefig(plotsdir("fields","reconstruction", "goal_marginalize_FRACTION_quantification_of_goal_missing_samples.svg"))
    Plots.heatmap(utils.squeeze(sum(Float64.(X.occzeroinds),dims=place_dims)), title="Quantification of missing samples\nin x and y of (x,y,γ,p)\n")
    Plots.savefig(plotsdir("fields","reconstruction", "place_marginalize_quantification_of_goal_missing_samples.svg"))
    Plots.heatmap(utils.squeeze(mean(Float64.(X.occzeroinds),dims=place_dims)), title="Quantification of missing samples\nin x and y of (x,y,γ,p)\n")
    Plots.savefig(plotsdir("fields","reconstruction", "place_marginalize_FRACTION_quantification_of_goal_missing_samples.svg"))

end
