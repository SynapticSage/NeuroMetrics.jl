quickactivate("/home/ryoung/Projects/goal-code/")
@time include(scriptsdir("fields", "Initialize.jl"))

# ----------------------------------------
# FUNCTION SPECIFIC SHORTCUTS AND SETTINGS
# ----------------------------------------
# Marginals and reconstructions, and data structures to map em out
shortcut_names = OrderedDict(
                             "currentAngle"=>"γ",
                             "currentPathLength"=>"p",
                             "stopWell"=>"G")
𝕀(d) = Dict(zip(values(shortcut_names), keys(shortcut_names)))[d]
reconstruction_comparisons = Dict( 
                                  "vs(angle,place)"            => ("γ|x,y","x,y|γ"),
                                  "vs(spect-angle,spec-place)" => ("γ,G|x,y","x,y|γ,G"),
                                  "vs(spec-angle,place)"       => ("p,γ,G|x,y,G","x,y,G|p,γ,G"),
                                  "vs(goal,place)"             => ("p,γ,G|x,y,G","x,y,G|p,γ,G"),
                                  "vs(spec-goal,place)"        => ("p,γ,G|x,y,G","x,y,G|p,γ,"),
                                  "vs(spec-goal,spec-place)"   => ("p,γ,G|x,y,G","x,y,G|p,γ,G"))
reconstructions_required = vec([x[i] for x in values(reconstruction_comparisons), i in 1:2])
props = ["x", "y", "currentPathLength", "currentAngle","stopWell"]

# Shortcut functions
"""
Returns the integer dim indices for a prop-string e.g. "x-y"->[1,2]
"""
function 𝔻(dimstr)
    out = [findfirst(dim.==dims) for dim in split(dimstr,",")]
    if out isa Vector{Nothing}
        @error "Nope! One your dims=$dimstr is wrong. Check for missing syntax (commas)"
    end
    return out
end
"""
Returns the remaining dimensions not covered by a prop-string
"""
𝔻̅(dimstr) = setdiff(1:length(props), 𝔻(dimstr)) # dims inverse
"""
Returns remaning dimensions as a joined prop string, instead of ints
"""
𝔻̅ⱼ(dimstr) = join(dims[𝔻̅(dimstr)],",") # joined
"""
Returns string with shortcut names
and the ₀ version returns the remaining names
"""
ℝ(dims) = replace(dims, [shortcut_names...][begin:1:end]...) # replace
ℝ₀ⱼ(dims) = ℝ(join(props[𝔻₀(dims)],"-")) # joined

dims  = ℝ(props)
x = Set(vec([split(x,"|")[1] for x in reconstructions_required])) # requires the LHS andd inverse of the RHS of each reconstruction
y = Set(vec([𝔻̅ⱼ(split(x,"|")[2]) for x in reconstructions_required])) # requires the LHS andd inverse of the RHS of each reconstruction
marginals_required = x ∪ y

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
utils.pushover("Processed up to joint distribution")

# ---------
# MARGINALS
# ---------
# Acquire marginals P(X,Y), P(γ, p, G)
@time @showprogress for marginal in marginals_required
    d̅ =  𝔻̅(marginal)
    println("marginal=>$marginal d̅ = $(d̅)")
    @time F[marginal] = operation.marginalize(X, dims = d̅ );
end
utils.pushover("Finished marginals")

# ---------------
# Reconstructions
# ---------------
# Obtain reconstructions!
R̂ = Dict()
@time for reconstruction in reconstructions_required
    given = split(reconstruction, "|")[2]
    inverse_given = join(dims[𝔻₀(given)], ",")
    @time R̂[reconstruction] = operation.apply(model.reconstruction, F["placegoal-joint"].occR, F[inverse_given].Rₕ);
end

# ---------------
# SUMMARIES
# ---------------
# Get reconstruction model error summary
E = Vector{DataFrame}([])
for reconstruction in reconstructions_required
    what, given = split(reconstruction, "|")
    error = model.reconstruction_error(F[what].Rₕsq, R̂[reconstruction])
    error = table.to_dataframe(error; name="error")
    append!(E, error)
end
E = vcat(error..., source=:model=>[reconstructions_required])

# Stacked summary!
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
