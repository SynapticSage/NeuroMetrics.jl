include("prep2plot.jl")

# ------------------------------------------------------------
# PLOTTING for DFS, time-averaged summary of DF
# ------------------------------------------------------------


dfs, dfsn = subset_dfs()
difference_in_react_match_nonmatch(dfs, dfsn)
# ------------------------------------------------
begin
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
        push!(h,plot(hh, title="Reactivation scores\nfor matched sets $areas\n\n\n\n", 
            xlabel="Reactivation score", ylabel="Count"))
    end
    plot(h..., layout=(1,3), size=(1200, 400))
    ylims!(0,2000)

    h=[]
    for epochs in unique(DF.epoch)
        dfs_match, dfs_nonmatch = subset_dfs(:epoch => x-> x .== epochs)
        hh=difference_in_react_match_nonmatch(dfs_match, dfs_nonmatch)
        push!(h,plot(hh, title="Reactivation scores for matched sets\nepoch=$epochs", 
                         xlabel="Reactivation score", ylabel="Count"))
    end
    plot(h..., layout=(1,3), size=(1200, 400))
end

begin
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

end

