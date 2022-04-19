quickactivate("/home/ryoung/Projects/goal-code/")
@time include(scriptsdir("fields", "Initialize.jl"))
using NaNStatistics
using StatsBase
using ColorSchemes

# ----------------------------------------
# FUNCTION SPECIFIC SHORTCUTS AND SETTINGS
# ----------------------------------------
# Marginals and reconstructions, and data structures to map em out
shortcut_names = OrderedDict(
                             "currentHeadEgoAngle"=>"γ",
                             "currentPathLength"=>"p",
                             "stopWell"=>"G")
si = operation.selectind
𝕀(d) = Dict(zip(values(shortcut_names), keys(shortcut_names)))[d]
recon_compare = Dict( 
      "vs(γ|ℙ,ℙ|γ)"     => ("γ|x,y","x,y|γ"),
      "vs(γₛ|ℙ,ℙ|γₛ)"   => ("γ,G|x,y","x,y|γ,G"),
      "vs(𝔾|ℙ,ℙ|𝔾)"     => ("p,γ|x,y","x,y|p,γ"),
      "vs(𝔾ₛ|ℙ,ℙ|𝔾ₛ)"   => ("p,γ,G|x,y","x,y|p,γ,G"),
      "vs(𝔾ₛ|ℙₛ,ℙₛ|𝔾ₛ)" => ("p,γ,G|x,y,G","x,y,G|p,γ,G"),
      "vs(ℙₛ|ℙ,ℙ|ℙₛ)"   => ("x,y,G|x,y","x,y|x,y,G"),
     )
inv(x) = Dict(zip(values(x), keys(x)))
recon_name(x, op) = replace(x, "vs("=>"ε(", ","=>") $op ε(")
recon_req = vec([x[i] for x in values(recon_compare), i in 1:2])
recon_req = [Set(recon_req)...]
props = ["x", "y", "currentPathLength", "currentHeadEgoAngle","stopWell"]

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
x = Set(vec([split(x,"|")[1] for x in recon_req])) # requires the LHS andd inverse of the RHS of each reconstruction
y = Set(vec([𝔻̅ⱼ(split(x,"|")[2]) for x in recon_req])) # requires the LHS andd inverse of the RHS of each reconstruction
marginals_required = x ∪ y

# ----------------
# Helper functions
# ----------------
function perform_reconstructions_marginals_and_error(beh, spikes, field_kws; 
        F::Dict=Dict())
    @time X = field.get_fields(beh, spikes; field_kws...);
    #F["placegoal-joint"] = X
    field_size = size(operation.selectind(X.Rₕ))

    # ---------
    # MARGINALS
    # ---------
    # Acquire marginals P(X,Y), P(γ, p, G)
    @time @showprogress for marginal in marginals_required
        d̅ =  𝔻̅(marginal)
        println("marginal=>$marginal d̅ = $(d̅)")
        @time F[marginal] = operation.marginalize(X, dims = d̅ );
    end

    # ---------------
    # Reconstructions
    # ---------------
    # Obtain reconstructions!
    R̂ = Dict()
    @time for reconstruction in recon_req
        dimr, given = split(reconstruction, "|")
        #inverse_given = join(dims[𝔻₀(given)], ",")
        #ig_set   = split(inverse_given,",")
        marginalize_dims = 𝔻̅(dimr)
        println(reconstruction)
        @time R̂[reconstruction] = operation.apply(model.reconstruction, 
                                                  X.occR, 
                                                  F[given].Rₕ;
                                                  marginalize_dims=marginalize_dims,
                                                 );
        @assert ndims(si(R̂[reconstruction])) == length(split(dimr,","))
        @assert all(size(si(R̂[reconstruction])) .== field_size[𝔻(dimr)])
    end

    # ---------------
    # SUMMARIES
    # ---------------
    # Get reconstruction model error summary
    E = Vector{DataFrame}([])
    for reconstruction in recon_req
        what, given = split(reconstruction, "|")
        error = model.reconstruction_error(F[what].Rₕsq, R̂[reconstruction])
        error = table.to_dataframe(error; name="error")
        push!(E, error)
    end
    E = vcat(E..., source=:model=>recon_req)
    what,under = [vec(x) for x in eachrow(cat(split.(E.model,"|")...; dims=2))]
    E.what, E.under = what, under
    E = sort(E, [:model, :area, :unit])
    utils.pushover("Finished reconstruction summaries")
    return E
end
function create_unstacked_error_table(E, recon_compare)
    uE = unstack(E[!,Not([:what,:under])], :model,:error)[!,Not([:dim_1, :dim_2])]
    uE.∑ε = vec(nansum(replace(Matrix(uE[:, recon_req]),missing=>NaN); dims=2))
    uE = sort(uE, [:area,:∑ε])
    for (rc, compare) in recon_compare
        #uE[!,rc*"div"] = uE[!, compare[1]] ./ uE[!, compare[2]]
        println(rc)
        uE[!,rc] = (uE[!, compare[1]] .- uE[!, compare[2]])./(uE[!,compare[1]] .+ uE[!,compare[2]])
        uE[!,rc] = (uE[!, compare[1]] .- uE[!, compare[2]])
    end
end

# ----------------------------------------
# PLACE-GOAL JOINT DISTRIBUTION P(x,y, γ,p,G)
# ----------------------------------------
filters = merge(kws.filters,
                filt.correct,
                filt.notnan("currentHeadEgoAngle"), 
                filt.minmax("currentPathLength", 2, 150))
field_kws = (; kws..., resolution=[40, 40, 40, 40, 5], gaussian=0, props=props,
          filters=merge(kws.filters, filters))
E = perform_reconstructions_marginals_and_error(beh, spikes, field_kws; F=F)

# ----------------------------------------
# PLACE-GOAL-headdir JOINT DISTRIBUTION P(x,y,H,G)
# ----------------------------------------
perform_reconstructions_marginals_and_error(beh, spikes, field_kws; F=F, E=E)

# Acquire tidy unstacked representation
uE = create_unstacked_error_table(E, recon_compare)
                    
#,---.|         |    
#|---'|    ,---.|--- 
#|    |    |   ||    
#`    `---'`---'`---'
if ploton

    function r(df)
        df = copy(df)
        reps = merge(Dict(n=>n for n in names(df) if !(occursin("vs",n))),
                     Dict(n=>replace(recon_name(n,"-")," "=>"") for n in names(df) if occursin("vs",n)))
        println(reps)
        rename(df, reps...)
    end
    _mod(df::DataFrame) = [x  for x in names(df) if occursin("|", x) && !occursin("ε", x) && !(occursin("vs",x))]
    _com(df::DataFrame) = [x for x in sort(names(df)) if (occursin("ε", x) || occursin("vs",x)) && !occursin("∑", x)]
    MOD(df::DataFrame) = df[!, _mod(df)]
    COM(df::DataFrame) = df[!, _com(df)]
    amM(df)  = (minimum(Matrix(df)), maximum(Matrix(df)))
    smM(df)  = (-maximum(abs.(Matrix(df))), maximum(abs.(Matrix(df))))
    ruE = r(uE)


    ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
    #           SUMMARY ..... TYPICAL ERROR VALUES ....    #
    ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
    histogram(vec(Matrix(uE[:,3:end])), bins=100, xticks=(0:0.2:0.8), yscale=:log10, label="typical errors")
    savefig(plotsdir("fields","reconstruction","hist_typical_errors.png"))


    ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
    #           SUMMARIES ... of reconstructions ...       #
    ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##


    function K(total, vert=true) 
        if vert
        return (;layout= grid(total,1), tickfontsize=3, legendfontsize=4, link=:x)
        else
        return (;tickfontsize=3, legendfontsize=4, link=:x)
        end
    end

    Δx = 0.00001
    # Create summary of reconstructions: Cumulative distributions
    P, AUC = [], []
    colors=get(colorschemes[:isoluminant_cgo_80_c38_n256], 1:length(_mod(uE)), :extrema)
    total = length(_mod(uE))
    for (i,col) in enumerate(_mod(uE))
        x = uE[:,col]
        x = ecdf(replace(x, missing=>NaN))
        vals=x.sorted_values[1]:Δx:x.sorted_values[end]
        color=colors[i]
        p = Plots.plot(vals, x(vals), fillrange=0, label=col, color=color)
        auc = cumsum(x(vals))[end]*Δx
        push!(AUC,auc)
        Plots.hline!([0.5], c=:black, linestyle=:dash, label="half mast")
        Plots.annotate!(quantile(vals,0.5), 0.5*1.2, text("auc=$(@sprintf("%2.2f",auc))", :black, 3))
        push!(P, p)
    end
    p1 = Plots.plot(P[sortperm(AUC)]...; K(total,false)..., margin=-1mm, legendposition=:bottomright)

    P, AUC = [], []
    total = length(_com(ruE))
    colors=get(colorschemes[:tableau_red_green_gold], 1:total, :extrema)
    for (i, col) in enumerate(_com(ruE))
        x = ruE[:,col]
        x = ecdf(replace(x, missing=>NaN))
        vals=x.sorted_values[1]:Δx:x.sorted_values[end]
        color=colors[i]
        auc = cumsum(x(vals))[end]*Δx.*sign.(x(vals))
        push!(AUC,auc)
        p = Plots.plot(vals, x(vals), fillrange=0, label=col, link=:x, c=color, color=color)
        Plots.vline!([0], c=:black, linestyle=:dash, label="no difference", legendposition=:left)
        #Plots.annotate!(quantile(vals,0.5), 0.5*1.2, text("auc=$(@sprintf("%2.2f",auc))", :black, 3))
        push!(P, p)
    end
    p2= Plots.plot(P[sortperm(AUC)]...; K(total, true)..., link=:x, margin=-1mm)

    Plots.plot(p1, p2)


    # Plot out errors and differences per cell
    layout = @layout [[a; b] c]
    Plots.plot(
    heatmap(_mod(ruE),1:size(ruE,1), Matrix(MOD(ruE)),  clims=amM(MOD(ruE)), xrotation=45),
    heatmap(_mod(ruE), _mod(ruE), cor(Matrix(MOD(ruE))), xrotation=45),
    heatmap(_com(ruE), 1:size(ruE,1), Matrix(COM(ruE)), clims=smM(COM(ruE)), xrotation=20, c=:vik),
    #heatmap(com(ruE), com(ruE), cor(Matrix(COM(ruE))), clims=(-1,1), xrotation=45, c=:vik),
    layout=layout)


    ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
    #              RECONSTRUCTIONS                         #
    ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##

    # 1D or 2D
    # Pick a model
    using Interact, Blink
    W = []
    for i in [1,2,4]
        C = _com(uE)[i]
        best = sort(uE[!, ["unit", "area", C]], C, rev=true)
        recon₁, recon₂ = recon_compare[C]
        marg₁, marg₂ = split(recon₁,"|")[1], split(recon₂,"|")[1]
        sk = operation.sk
        ski(f::Dict, p::Pair...) = si(sk(f, p...))
        ui = @manipulate for neuron in best.unit
            println("neuron=$neuron")
            f₁, f₂ = ski(F[marg₁].Rₕsq, "unit"=>neuron), ski(F[marg₂].Rₕsq, "unit"=>neuron)
            r₁, r₂ = ski(R̂[recon₁], "unit"=>neuron), ski(R̂[recon₂], "unit"=>neuron)
            Plots.plot(
                       field.plot.show_field(f₁, key=(;n=neuron, m=marg₁), textcolor=:black, fontsize=10),
                       field.plot.show_field(f₂, key=(;n=neuron, m=marg₂), textcolor=:black, fontsize=10),
                       field.plot.show_field(r₁, key=(;n=neuron, m=recon₁),textcolor=:black, fontsize=10),
                       field.plot.show_field(r₂, key=(;n=neuron, m=recon₂),textcolor=:black, fontsize=10),
                      )
        end
        w = Window()
        body!(w,ui)
        push!(W,w)
    end


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


    
    ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
    # Applet for visualizing top difference comparisons    #
    ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##



    ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
    #                   SUMMARY OF MISSING SAMPLES         #
    ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
    goal_dims = [3,4,5]
    place_dims = [1,2]
    Plots.heatmap(utils.squeeze(sum(Float64.(X.occzeroinds),dims=goal_dims)),
                  title="Quantification of missing samples\nin γ and p of (x,y,γ,p)\n")
    Plots.savefig(plotsdir("fields","reconstruction", "goal_marginalize_quantification_of_goal_missing_samples.svg"))
    Plots.heatmap(utils.squeeze(mean(Float64.(X.occzeroinds),dims=goal_dims)),
                  title="Quantification of missing samples\nin γ and p of (x,y,γ,p)\n")
    Plots.savefig(plotsdir("fields","reconstruction", "goal_marginalize_FRACTION_quantification_of_goal_missing_samples.svg"))
    Plots.heatmap(utils.squeeze(sum(Float64.(X.occzeroinds),dims=place_dims)),
                  title="Quantification of missing samples\nin x and y of (x,y,γ,p)\n")
    Plots.savefig(plotsdir("fields","reconstruction", "place_marginalize_quantification_of_goal_missing_samples.svg"))
    Plots.heatmap(utils.squeeze(mean(Float64.(X.occzeroinds),dims=place_dims)),
                  title="Quantification of missing samples\nin x and y of (x,y,γ,p)\n")
    Plots.savefig(plotsdir("fields","reconstruction", "place_marginalize_FRACTION_quantification_of_goal_missing_samples.svg"))

end


# ADD DATA TO CELL STRUCTURE FOR OTHER ANALYSES
cells = leftjoin(cells, ue[:,[:unit,]])
