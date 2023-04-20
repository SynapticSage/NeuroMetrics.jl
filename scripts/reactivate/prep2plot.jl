using DrWatson
DrWatson.quickactivate("/home/ryoung/Projects/goal-code-run/")
cd(scriptsdir("reactivate"))
DrWatson.quickactivate("/home/ryoung/Projects/goal-code")
include("imports.jl") # include(scriptsdir("reactivate", "imports.jl"))

if !isdefined(Main, :opt)
    include("run.jl") 
end



# ------------------------------------------------
#               HELPER FUNCTIONS
# ------------------------------------------------
"""
get_pmatch_cols(decouple_cols)
Return the columns to partially match the actual labels of the samples (what
animals actually did) against the labels of the templates (the reactivation
template that was used).
# Globals used
- `default_match_cols`: Default columns to match.
# Arguments
- `decouple_cols`: Columns to decouple from the default matching columns.
"""
function get_pmatch_cols(decouple_cols=[])
    [
    [x for x in default_match_cols[1] if 
        all(.!occursin.(string.(decouple_cols), string(x)))
    ],
    [x for x in default_match_cols[2] if 
        all(.!occursin.(string.(decouple_cols), string(x)))
        ]
    ]
end
"""
set_pmatch!(decouple_cols)
Set the `:pmatch` column in `DFS` to `true` if 
are equal.
# Globals used
- `DFS`: Dataframe to set the `:pmatch` column in.
- `default_match_cols`: Default columns to match.
# Arguments
- `decouple_cols`: Columns to decouple from the default matching columns.
"""
function set_pmatch!(decouple_cols)
    throw(Error("set_pmatch! causes huge memory usage..."))
    nothing
end
"""
    subset_dfs(pos...)
Subset the dataframe `DFS` by the conditions in `general_conditions` and
`pos...`. If `:exclude` is defined in the subset, then exclude those rows
where `:exclude` is `true`.
We generate two subsets of `DFS` where `:match` is `true` and `false`,
respectively. Match 
# Globals used
- `DFS`: Dataframe to subset.
- `general_conditions`: Conditions to subset `DFS` by.
# Arguments
- `pos...`: Conditions to subset `DFS` by.
"""
function subset_dfs(pos...)
    # Subset of DF where k_tmpl == k_test
    dfs_match    = subset(DFS, :pmatch => x -> x .== true, 
                          enforce_conditions..., pos...)
    dfs_nonmatch = subset(DFS, :pmatch => x -> x .== false,
                          enforce_conditions..., pos...)
    # If exclude defined in column
    if :exclude in propertynames(dfs_match)
        dfs_match =    subset(dfs_match,    :exclude => x -> x .== false)
        dfs_nonmatch = subset(dfs_nonmatch, :exclude => x -> x .== false)
    end
    return dfs_match, dfs_nonmatch
end
# ------------------------------------------------
function difference_in_react_match_nonmatch(dfs_match, dfs_nonmatch)
    # do matching sets have a higher reactivation score?
    begin
        println("Matched sets have a higher reactivation score than non-matched sets?")
        # HypothesisTests.UnequalVarianceTTest(df_match.value, 
        #                                      df_nonmatch.value)
        println("DFS:")
        HypothesisTests.UnequalVarianceTTest(dfs_match.mean, 
                                             dfs_nonmatch.mean)
    end
    h=histogram(dfs_match.mean, label="Matched", 
                linewidth=0, strokealpha=0.01,
              title="Reactivation scores for matched sets", alpha=0.5)
    histogram!(dfs_nonmatch.mean, label="Non-matched", 
                linewidth=0, strokealpha=0.01,
               title="Reactivation scores for non-matched sets", alpha=0.5)
    vline!([0], label="", linecolor=:black, linestyle=:dash)
    vline!([mean(dfs_match.mean)], label="Mean matched", 
           linecolor=:blue, linestyle=:dashdot)
    vline!([mean(dfs_nonmatch.mean)], label="Mean non-matched",
            linecolor=:red, linestyle=:dashdot)
    xlims!(-0.2, 1)
    ylims!(0, 10000)
    return h
end
"""
    match_nonmatch_splits(split_tmpl_match = [:startWell],
                          template_dims = [:startWell_tmpl, :stopWell_tmpl])
Measure the difference in reactivation scores between tmpl and test
sets (we iterate these). Additionally, the .pmatch property of the dataframe
operated upon is used to further split these data into two result pools:
"matched" and "non-matched".
I may need to play around with the semantics of these. It's confusing the hear
the word "match" in two contexts within this function.
# Arguments
- `split_actual_dims`: which actual dimensions to split on (ie what actually
happened to the animal when measurements are taken)
- `template_dims`: which template dimensions to split on (ie which reactivation
template was used to measure against)
# Returns
- `matching`: an array of dataframe cuts where .pmatch column is true
- `nonmatching`: an array of dataframe cuts where .pmatch column is false
- `tml_matches` : whether the `actual_indices` match the `template_indices`
                  per entry
"""
function match_nonmatch_splits(actual_dims = [:startWell],
        template_dims = [:startWell_tmpl, :stopWell_tmpl])
    # Step 0 : Handle input types
    actual_dims = actual_dims isa Symbol ? [actual_dims] : actual_dims
    template_dims = template_dims isa Symbol ? [template_dims] : template_dims
    # Step 1 : Match all properties except unmatch variable
    dfs, _ = subset_dfs()
    sort!(dfs, [template_dims..., actual_dims...])
    U = Dict(x=>unique(dfs[!, x]) for x in propertynames(dfs))
    #ISSUE: why U[:startWell] has only 4 values?
    match_dims = Symbol.(replace.(string.(template_dims),"_tmpl"=>"")) .==
                 (actual_dims...,)
    # Step 2 : Prepare our loop variables
    tmpl_iters = Iterators.product([U[x] for x in template_dims]...)
    tmpl = first(tmpl_iters)
    actual_iters = Iterators.product([U[x] for x in actual_dims]...)
    actual_iter_indices = collect(Iterators.product(axes(actual_iters)...))
    matching    = Array{Union{DataFrame,Missing}}(missing,
    size(tmpl_iters)..., size(actual_iters)...)
    nonmatching = Array{Union{DataFrame, Missing}}(missing,
    size(tmpl_iters)..., size(actual_iters)...)
    tmpl_matches = Array{Union{Missing,Bool}}(missing, size(tmpl_iters)..., 
        size(actual_iters)...)
    tmp_iter_indices = collect(Iterators.product(axes(tmpl_iters)...))
    for (tmpl_ind, tmpl) in zip(tmp_iter_indices, tmpl_iters)
        selector_tmpl = [template_dims[i] => 
                            x -> x .== tmpl[i] for i in 1:length(tmpl)]
        dfsm, dfsn = subset_dfs( selector_tmpl...)
        sort!(dfsm, [template_dims..., actual_dims...])
        sort!(dfsn, [template_dims..., actual_dims...])
        (i, actual) = (1,first(actual_iters))
        for (actual_ind, actual) in zip(actual_iter_indices, actual_iters)
            selector_te = [actual_dims[i] => x -> x .== actual[i] for i in
                            1:length(actual)]
            dm = subset(dfs_match, selector_te...)
            dn = subset(dfs_nonmatch, selector_te...)
            matching[tmpl_ind..., actual_ind...]    = dm
            nonmatching[tmpl_ind..., actual_ind...] = dn
            tmpl_matches[tmpl_ind..., actual_ind...] = 
                        all(tmpl[match_dims] .== actual)
        end
    end
    matching = DimArray(matching, (template_dims..., actual_dims...))
    nonmatching = DimArray(nonmatching, (template_dims..., actual_dims...))
    return matching, nonmatching, tmpl_matches
end
"""
    tmpl_ylabels(gdf::GroupedDataFrame)   
get the ylabels for a grouped dataframe based on the template
"""
function tmpl_ylabels(gdf::GroupedDataFrame)
    map(tmpl_ylabels, gdf |> collect)
end
"""
    get_tmpl_ylabels(df::DataFrame)
get the ylabels for a dataframe based on the template
"""
function tmpl_ylabels(df::Union{SubDataFrame, AbstractDataFrame})
    df=first(eachrow(df))
    actual = "$(df.startWell)-$(df.stopWell)"
    label  = "$(df.startWell_tmpl)-$(df.stopWell_tmpl)"
    if df.i_tmpl == 0 
        "OTHER"
    else
        actual == label ? "ACTUAL: $label" : label
    end
end
"""
    tmpl_labels_dict(df::DataFrame, type="tmpl")
Return a dictionary of labels for the template
"""
function tmpl_labels_dict(df::DataFrame, type="tmpl"; move_state=false)
    index     = type == "tmpl" ? :i_tmpl : :i_test
    startWell = type == "tmpl" ? :startWell_tmpl : :startWell
    stopWell  = type == "tmpl" ? :stopWell_tmpl : :stopWell
    moving    = type == "tmpl" ? :moving_tmpl : :moving
    labels = Dict()
    df = vcat(DataFrame.(
        unique(eachrow(df[!,[index, startWell, stopWell, moving]]))
    )...)
    # replace!(df[!,moving], Dict(true=>"moving", false=>"still"))
    for (i, s, S, ms) in zip((df[!, index]), (df[!, startWell]),
                    (df[!, stopWell]), (df[!, moving]))
        ms = ms == true ? "M" : "S"
        if i != 0
            labels[i] = move_state === false ? "$s-$S" : "$ms: $s-$S"
        else
            labels[i] = "OTHER"
        end
    end
    return labels
end
"""
    plot_tmpl_match(df::DataFrame)
Return a fillstyle for a given value of correct
"""
function corerr_fillstyle(iscorrect)
    if iscorrect == 1
        return nothing
    elseif iscorrect == 0
        return :\
    else
        return :x
    end
end


# DF = DIutils.arr.get_quantile_filtered(DF, :value, 0.001)
# Want to select train=moving and test=stationary
enforce_conditions = [:moving_tmpl =>  x -> x .== true,
                      :moving       => x -> x .== false]
default_match_cols =  [[:startWell, :stopWell, :ha, :epoch],
                       [:startWell_tmpl, :stopWell_tmpl, :ha_tmpl, :epoch_tmpl]]
match_cols = get_pmatch_cols()
println("Matching columns: ", match_cols)
DFS[!,:pmatch] .= all(Matrix(DFS[!, match_cols[1]]) .== Matrix(DFS[!, match_cols[2]]), dims=2) |> vec
DFS.exclude = DFS.n .< 10

if !isdefined(Main, :beh)
    beh = DI.load_behavior(opt["animal"], opt["day"])
end
