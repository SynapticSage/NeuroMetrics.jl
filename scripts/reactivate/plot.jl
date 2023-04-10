include("run.jl") # Just so language-server protocol can find the symbols
                  # from the script that generated the data.

# ------------------------------------------------------------
#  _| || |_  | | | |_   _ _ __   | |_ ___  ___| |_(_)_ __   __ _ 
# |_  ..  _| | |_| | | | | '_ \  | __/ _ \/ __| __| | '_ \ / _` |
# |_      _| |  _  | |_| | |_) | | ||  __/\__ \ |_| | | | | (_| |
#   |_||_|   |_| |_|\__, | .__/   \__\___||___/\__|_|_| |_|\__, |
#                   |___/|_|                               |___/ 
# ------------------------------------------------------------
DF = DIutils.arr.get_quantile_filtered(DF, :value, 0.001)
# Want to select train=moving and test=stationary
general_conditions = [:moving_tmpl => x -> x .== true,
                      :moving       => x -> x .== false]

default_match_cols =  [[:startWell, :stopWell, :startstopWell, :ha, :epoch],
[:startWell_tmpl, :stopWell_tmpl, :startstopWell_tmpl, :ha_tmpl, :epoch_tmpl]]
match_cols = copy(default_match_cols)
DFS.match = all(Matrix(DFS[!, match_cols[1]]) .== Matrix(DFS[!, match_cols[2]]),
                dims=2) |> vec
DFS.exclude = DFS.n .< 10


"""
    subset_dfs(pos...)
Subset the dataframe `DFS` by the conditions in `general_conditions` and
`pos...`. If `:exclude` is defined in the subset, then exclude those rows
where `:exclude` is `true`.
We generate two subsets of `DFS` where `:match` is `true` and `false`,
respectively. Match 
# Arguments
- `pos...`: Conditions to subset `DFS` by.
"""
function subset_dfs(pos...)
    # Subset of DF where k_tmpl == k_test
    dfs_match    = subset(DFS, :match => x -> x .== true, 
                          general_conditions..., pos...)
    dfs_nonmatch = subset(DFS, :match => x -> x .== false,
                          general_conditions..., pos...)
    # If exclude defined in column
    if :exclude in propertynames(dfs_match)
        dfs_match =    subset(dfs_match,   :exclude => x -> x .== false)
        dfs_nonmatch = subset(dfs_nonmatch, :exclude => x -> x .== false)
    end
    return dfs_match, dfs_nonmatch
end
dfs_match, dfs_nonmatch = subset_dfs()

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
difference_in_react_match_nonmatch(dfs_match, dfs_nonmatch)
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
# Step 1: Change .match to encode match excepting the startWell
# Step 2: function that subsets, hypothesis tests to the REPL, and
#         plots the histogram as above
# ------------------------------------------------
# Step 1
unmatch = :startWell
match_cols = [[x for x in default_match_cols[1] if 
                !occursin(string(unmatch), string(x))],
              [x for x in default_match_cols[2] if 
                !occursin(string(unmatch), string(x))]]
DFS.match = all(Matrix(DFS[!, match_cols[1]]) .== Matrix(DFS[!, match_cols[2]]),
                dims=2) |> vec;
DFS.exclude = DFS.startWell == -1 || DFS.stopWell == -1 || 
              DFS.startWell_tmpl == -1 || DFS.stopWell_tmpl == -1
dfs, _ = subset_dfs()
U = Dict(x=>unique(dfs[!, x]) for x in propertynames(dfs))
train_dims = [:startWell_tmpl, :stopWell_tmpl]
test_dims = [:startWell]
# Step 2
h = []
train_iters = Iterators.product([U[x] for x in train_dims]...)
item = first(train_iters)
for item in train_iters
    selector_tr = [train_dims[i] => x -> x .== item[i] for i in 1:length(item)]
    dfs_match, dfs_nonmatch = subset_dfs( selector_tr...)
    test_iters = Iterators.product([U[x] for x in test_dims]...)
    matching    = Array{Any}(undef, size(train_iters)..., size(test_iters)...)
    nonmatching = Array{Any}(undef, size(train_iters))
    startWell = first(test_iters)
    for (i,startWell) in enumerate(unique(dfs.startWell))
        selector_te = [test_dims[i] => x -> x .== startWell for i in 1:length(startWell)]
        dm = subset(dfs_match, selector_te...)
        dn = subset(dfs_nonmatch, selector_te...)
        if startWell_tmpl == startWell
            push!(matching, dm)
        else
            push!(nonmatching, dn)
        end
    end
end




