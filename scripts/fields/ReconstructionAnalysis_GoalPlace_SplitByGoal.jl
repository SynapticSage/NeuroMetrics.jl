quickactivate("/home/ryoung/Projects/goal-code/")
@time include(scriptsdir("fields", "Initialize.jl"))
recon = field.recon
#recon_process =  field.recon_process
using NaNStatistics
using StatsBase
using ColorSchemes
using Printf
using ThreadSafeDicts

F = ThreadSafeDict(pairs(F)...)

# ---------------------------------------- FUNCTION SPECIFIC SHORTCUTS AND SETTINGS ----------------------------------------
si,sk = operation.selectind, operation.selectkey
filters = merge(kws.filters,
                filt.correct,
                filt.notnan("currentHeadEgoAngle"), 
                filt.minmax("currentPathLength", 2, 150))
overall_compare = Dict()


# ----------------------------------------
# PLACE-GOAL JOINT DISTRIBUTION P(x,y, γ,p,G)
# ----------------------------------------
recon_compare = Dict( 
      "vs(γ|ℙ,ℙ|γ)"     => ("γ|x,y","x,y|γ"),
      "vs(γₛ|ℙ,ℙ|γₛ)"   => ("γ,G|x,y","x,y|γ,G"),
      "vs(𝔾|ℙ,ℙ|𝔾)"     => ("p,γ|x,y","x,y|p,γ"),
      "vs(𝔾ₛ|ℙ,ℙ|𝔾ₛ)"   => ("p,γ,G|x,y","x,y|p,γ,G"),
      "vs(𝔾ₛ|ℙₛ,ℙₛ|𝔾ₛ)" => ("p,γ,G|x,y,G","x,y,G|p,γ,G"),
      "vs(ℙₛ|ℙ,ℙ|ℙₛ)"   => ("x,y,G|x,y","x,y|x,y,G"),
     )
overall_compare = merge(overall_compare, recon_compare)
field_kws = (; kws..., resolution=[40, 40, 40, 40, 5], gaussian=0, 
             props=["x", "y", "currentPathLength",
                    "currentHeadEgoAngle","stopWell"],
          filters=merge(kws.filters, filters))
E = recon_process.perform_reconstructions_marginals_and_error(beh, spikes, field_kws; 
                                                              F=F,
                                               recon_compare=recon_compare)

# ----------------------------------------
# PLACE-GOAL-headdir JOINT DISTRIBUTION P(x,y,H,G)
# ----------------------------------------
recon_compare = Dict("vs(ℙₕ|ℙₛ,ℙₛ|ℙₕ)" => ("x,y,G|x,y,H","x,y,H|x,y,G"))
overall_compare = merge(overall_compare, recon_compare)
headdir_kws = (;field_kws..., props = ["x", "y", "headdir","stopWell"])
headdir_kws = (;headdir_kws..., resolution = [40,40,5,5])
E = recon_process.perform_reconstructions_marginals_and_error(beh, spikes, headdir_kws; 
                                                F=F, recon_summary=E, recon_compare);

# Acquire tidy unstacked representation
uE = recon_process.create_unstacked_error_table(E, recon_compare)

                    
#,---.|         |    
#|---'|    ,---.|--- 
#|    |    |   ||    
#`    `---'`---'`---'
if ploton

    function r(df)
        df = copy(df)
        reps = merge(Dict(n=>n for n in names(df) if !(occursin("vs",n))),
                     Dict(n=>replace(recon_process.get_recon_name(n,"-")," "=>"") for n in names(df) if occursin("vs",n)))
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
