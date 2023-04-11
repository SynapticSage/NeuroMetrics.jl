using DrWatson
DrWatson.quickactivate("/home/ryoung/Projects/goal-code-run/")
cd(scriptsdir("reactivate"))
DrWatson.quickactivate("/home/ryoung/Projects/goal-code")
include("imports.jl") # include(scriptsdir("reactivate", "imports.jl"))
include("run.jl") 
# Just so language-server protocol can find the symbols
# from the script that generated the data.
load_react_vars()
using DimensionalData


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
              title="Reactivation scores for matched sets", alpha=0.5)
    histogram!(dfs_nonmatch.mean, label="Non-matched", 
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
    tmpl_matches = Array{Union{Missing,Bool}}(missing, size(tmpl_iters)..., size(actual_iters)...)
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


# ------------------------------------------------------------
#  _| || |_  | | | |_   _ _ __   | |_ ___  ___| |_(_)_ __   __ _ 
# |_  ..  _| | |_| | | | | '_ \  | __/ _ \/ __| __| | '_ \ / _` |
# |_      _| |  _  | |_| | |_) | | ||  __/\__ \ |_| | | | | (_| |
#   |_||_|   |_| |_|\__, | .__/   \__\___||___/\__|_|_| |_|\__, |
#                   |___/|_|                               |___/ 
# ------------------------------------------------------------
DF = DIutils.arr.get_quantile_filtered(DF, :value, 0.001)
# Want to select train=moving and test=stationary
enforce_conditions = [:moving_tmpl =>  x -> x .== true,
                      :moving       => x -> x .== false]
default_match_cols =  [[:startWell, :stopWell, :ha, :epoch],
                       [:startWell_tmpl, :stopWell_tmpl, :ha_tmpl, :epoch_tmpl]]
match_cols = get_pmatch_cols()
println("Matching columns: ", match_cols)
DFS[!,:pmatch] .= all(Matrix(DFS[!, match_cols[1]]) .== Matrix(DFS[!, match_cols[2]]), dims=2) |> vec
DFS.exclude = DFS.n .< 10

dfs, dfsn = subset_dfs()
difference_in_react_match_nonmatch(dfs, dfsn)
# ------------------------------------------------

h = []
for k in 1:8
    dfs_match, dfs_nonmatch = subset_dfs(:component => x-> x .== k)
    hh = difference_in_react_match_nonmatch(dfs_match, dfs_nonmatch)
    push!(h, 
    plot(hh, xlabel="Reactivation score", ylabel="Count"))
end
plot(h..., layout=(2,4), size=(1200, 800), 
title="Reactivation scores for component 1-8")


h=[]
for areas in ["ca1-pfc", "ca1-ca1", "pfc-pfc"]
    dfs_match, dfs_nonmatch = subset_dfs(:areas => x-> x .== areas)
    hh=difference_in_react_match_nonmatch(dfs_match, dfs_nonmatch)
    push!(h,plot(hh, title="Reactivation scores for matched sets $areas", 
        xlabel="Reactivation score", ylabel="Count"))
end
plot(h..., layout=(1,3), size=(1200, 400))

h=[]
for epochs in unique(DF.epoch)
    dfs_match, dfs_nonmatch = subset_dfs(:epoch => x-> x .== epochs)
    hh=difference_in_react_match_nonmatch(dfs_match, dfs_nonmatch)
    push!(h,plot(hh, title="Reactivation scores for matched sets\nepoch=$epochs", 
                     xlabel="Reactivation score", ylabel="Count"))
end
plot(h..., layout=(1,3), size=(1200, 400))

# ------------------------------------------------
# Same as above section, but split out by different startWell to the
# same stopWell given (train_startWell, train_stopWell)
#
# In other words, we want to collect FOR each matching set (matched 
# on stopWell, ha, epoch) the reactivation score for each startWell
#
# Step 1: Change .pmatch to encode match excepting the startWell
# Step 2: function that subsets, hypothesis tests to the REPL, and
#         plots the histogram as above
# ------------------------------------------------
DFS.exclude = DFS.startWell      .== -1 .|| DFS.stopWell .== -1 .|| 
              DFS.startWell_tmpl .== -1 .|| DFS.stopWell_tmpl .== -1
decouple_cols = [:startWell]
match_cols = get_pmatch_cols(decouple_cols)
@assert(!(decouple_cols ⊆ match_cols[1]), "decouple_cols must not be in match_cols")
println("Matching columns: ", match_cols)
DFS[!,:pmatch] = all(Matrix(DFS[!, match_cols[1]]) .== Matrix(DFS[!, match_cols[2]]),
                dims=2) |> vec
println("Fraction of matches: ", mean(DFS.pmatch))
dfs,dfsn = subset_dfs()

# ISSUE: houston, we have a problem.  why is startwell_tmpl always ==
#         startwell?  
#------------------------------------------------
@assert .!all(DFS.startWell_tmpl .== DFS.startWell)
@assert .!all(dfs.startWell_tmpl .== dfs.startWell)
@assert .!all(dfsn.startWell_tmpl .== dfsn.startWell)
dfss = combine(groupby(dfs, [:startWell_tmpl, :stopWell_tmpl, :startWell]),
    :mean => mean, :mean => std, :mean => length)
sort!(dfss, [:stopWell_tmpl, :startWell_tmpl, :startWell, ])
dfss
#ISSUE:----------------------------------------------------------------

# Measure for template=[startWell, stopWell] and actual=[startWell]
@time matching, nonmatching, tmpl_matches = match_nonmatch_splits()
print("tmpl_matches: ");
tmpl_matches

# Now we want to summarize the effects and plot them out
tmpl_matches
μM = DimArray(map(collect(matching)) do x
    if x === missing
        missing
    else
        mean(x.mean)
    end
end, matching.dims)
μM = map(collect(matching)) do x
    if x === missing
        missing
    else
        mean(x.startWell)
    end
end
μN = map(collect(nonmatching)) do x
    if x === missing
        missing
    else
        mean(x.mean)
    end
end

heatmap.(eachslice(μN, dims=3))
#: ISSUE: THESE LOOK WRONG

# TRYING A SEPARATE APPROACH
# ------------------------------------------------
GC.gc()
match_cols = get_pmatch_cols([:startWell,:stopWell])
DF[!,:pmatch] = all(Matrix(DF[!, match_cols[1]]) .== Matrix(DF[!, match_cols[2]]),
                dims=2) |> vec
DFc = @subset(DF, :areas .== "ca1-ca1", :ha .== 'A', 
              :component .<= 5)
GC.gc()

using GoalFetchAnalysis.Munge.tensor
T = tensor.tensorize(DFS, [:startWell, :stopWell, :startWell_tmpl, :pmatch, :areas], [:mean])
heatmap.(eachslice(T[:, :, :, 1],dims=3))

#: TODO: RUN THIS
TT = nothing
GC.gc()
TT = tensor.tensorize(DFc, [:traj, :startWell, :stopWell, :startWell_tmpl, :stopWell_tmpl], [:time,:value])
#: BUG: WHY ARE THERE ONLY 12 trajectories?
