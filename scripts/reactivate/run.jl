if defined(Main, :DrWatson)
    include(scriptdir("reactivation","imports.jl"))
else
    include("imports.jl")
end

# Get and prep the data!
# ------------------------------------------------------------------
    spikes = DI.load_spikes(opt["animal"],   opt["day"])
    beh    = DI.load_behavior(opt["animal"], opt["day"])
    cells  = DI.load_cells(opt["animal"],    opt["day"])
    # Get firing rate DimArray
    R = spiking.torate(spikes, beh, gaussian=0.20)
    # Create dataframe of R
    Rdf = DataFrame(R, name=:rate)
    DIutils.pushover("Finished creating Rdf")
    # Create a new field enumerating (startWell, stopWell) combinations
    # with integers
    S  = eachrow(Matrix(beh[:, [:startWell, :stopWell]])) |> collect
    uS = unique(S, dims=1)
    sf = [findfirst(isequal(s), uS) for s in S]
    beh[:, :startstopWell] = sf
    # Create a field specifying whether the animal is still or moving
    beh[:, :moving] = beh[:, :speed] .> 2 # cm/s
    groupings = [:startstopWell, :startWell, :stopWell, :moving, :ha, :correct, 
        :epoch]
    # Register that to Rdf
    register(beh, Rdf,   on="time", transfer=["traj",string.(groupings)...])
    Rdf = @subset(Rdf,  
                        :stopWell .!= -1,
                        :startWell .!= 1)
    register(cells, Rdf, on="unit", transfer=["area"])
    ca1, pfc = groupby(Rdf, :area)[1:2];
    @assert(ca1.area[1] == "CA1" && pfc.area[1] == "PFC", 
    "Areas are not CA1 and PFC")
    ca1 = ca1[:, Not(:area)];
    pfc = pfc[:, Not(:area)];
    Rdf = Rdf[:, Not(:area)];
    dropmissing!(ca1, :startstopWell)
    dropmissing!(pfc, :startstopWell)
    DIutils.pushover("Finished creating area splits")
    ca1.ha = replace(ca1.ha, missing => 'm')
    ca1    = groupby(ca1, groupings);
    pfc    = groupby(pfc, groupings);
# ------------------------------------------------------------------
# Which data will we accept?, right now I'm only accepting combinations
# that have more than num_cells samples
# ------------------------------------------------------------------
# BUG:
# why so many 'h' 'a' missing? especially when start and stopwell are both
# defined? 
# It turns out ALL of the missing values hail from the incorrect trials
# ------------------------------------------------------------------
accepted = begin
    accepted = 
        Base.map(Base.filter(map(x->groupby(x,:unit)|>first,collect(ca1))) do x
            size(x,1) > size(cells,1) && length(unique(x.traj)) ≥ 3
    end) do x
        sz = size(x,1)
        ntraj = length(unique(x.traj))
        x=DataFrame(x[1,groupings])
        cols=propertynames(x)
        x[!,:size] .= sz
        x[!,:ntraj] .= ntraj
        x[!, [:size, :ntraj, cols...]]
    end
    accepted = sort(vcat(accepted...), :size, rev=true)
    accepted = dropmissing(accepted)
    # accepted |> Voyager()
        p=@df accepted scatter(:correct, :ha)
        p.subplots[1].series_list[1][:x] .+= 
                    0.02randn(size(p.subplots[1].series_list[1][:x]))
        p.subplots[1].series_list[1][:y] .+= 
                    0.05randn(size(p.subplots[1].series_list[1][:y]))
        plot(p)
    unicodeplots()
        histogram((accepted.size * 0.033333)./60)
        histogram(accepted.ntraj, xlims=(0,20))
        histogram(combine(groupby(accepted,[:startstopWell,:epoch]),
                :ntraj=>sum).ntraj_sum)
    gr()
    accepted
end

# ------------------------------------------
# Recombine ca1 and pfc, filter, and re-split
# ------------------------------------------
@time begin
    ca1 = combine(ca1, identity);
    pfc = combine(pfc, identity);
    inds = collect(eachrow(Matrix(ca1[!,groupings]))) .∈ 
          (collect(eachrow(Matrix(accepted[!,groupings]))),);
    inds = @time disallowmissing(replace(inds, missing=>false));
    print("Fraction of points kept: ", sum(inds)/length(inds))
    ca1  = ca1[inds,:];
    inds = collect(eachrow(Matrix(pfc[!,groupings]))) .∈ 
          (collect(eachrow(Matrix(accepted[!,groupings]))),);
    inds = disallowmissing(replace(inds, missing=>false))
    print("Fraction of points kept: ", sum(inds)/length(inds))
    pfc = pfc[inds,:];
    ca1 = groupby(ca1, groupings);
    pfc = groupby(pfc, groupings);
    if size(ca1,1) != size(pfc,1)
        println("size(ca1,1)=$(size(ca1,1)) != size(pfc,1)= $(size(pfc,1))")
        println("...Aligning!")
        K1, K2 = keys(ca1), keys(pfc)
        K = intersect(K1, K2)
        println("...Keeping $(length(K)) keys")
        println("---------------")
        println(K)
        println("---------------")
        pfc′, ca1′ = Dict(), Dict()
        for k in K
            k = NamedTuple(k)
            pfc′[k] = pfc[k]
            ca1′[k] = ca1[k]
        end
        ca1 = vcat(values(ca1′)...)
        pfc = vcat(values(pfc′)...)
    end
    ca1, pfc = groupby(ca1, groupings), groupby(pfc, groupings);
    @assert(size(ca1,1) == size(pfc,1),
    "size(ca1,1)=$(size(ca1,1)) != size(pfc,1)=$(size(pfc,1))")
    # @assert(map(
    # (x,y) -> size(x,1) == size(y,1), collect(ca1), collect(pfc)),
    # )
end

# ----------------------------------
# Now convert these back to DimArrays
# ----------------------------------
ca1 = map(zip(keys(ca1),ca1) |> collect) do (k,x)
    NamedTuple(k) => unstack(x, [:time,groupings...], 
        :unit, :rate; combine=mean)
end
pfc = map(zip(keys(pfc),pfc) |> collect) do (k,x)
    NamedTuple(k) => unstack(x, [:time,groupings...], 
        :unit, :rate; combine=mean)
end
ca1, pfc = Dict(ca1), Dict(pfc)

# ---------------------------
# Get props of each area split
# ---------------------------
startstop = map(values(ca1) |> collect) do x
    x.startstopWell[1]
end
startstop_pfc = map(values(pfc) |> collect) do x
    x.startstopWell[1]
end
@assert startstop == startstop_pfc
props = [groupings..., :time]
props_ca1 = map(ca1 |> collect) do (key,x)
    key=>NamedTuple([k=>x[!,k] for k in props])
end
props_pfc = map(pfc |> collect) do (key,x)
    key=>NamedTuple([k => x[!,k] for k in props])
end
props_ca1, props_pfc = Dict(props_ca1), Dict(props_pfc)

# ----------------------------
# Now we can test reactivation
# ----------------------------
Z_ca1, Z_pfc = Dict(), Dict()
for k in keys(ca1)
    Z_ca1[k] = reactivation.zscoreFRmatrix(ca1[k])
    Z_pfc[k] = reactivation.zscoreFRmatrix(pfc[k])
end

# ----------------------------
# Testing: DONE: PCAICA
# ----------------------------
begin
    zca1, zpfc = Z_ca1[1], Z_pfc[1]
    m = M_PCAICA()
    train = reactivation.TrainReact(m, zca1, zca1)
    m = ingredients(m, train)
    r1 = reactscore(m, zca1, zca1)
    m = M_PCAICA()
    train = reactivation.TrainReact(m, zca1, zpfc)
    m = ingredients(m, train)
    r2 = reactscore(m, zca1, zpfc)

    heatmap(replace(cor(r1), NaN=>0), c=:vik, clim=(-1,1))
    heatmap(replace(cor(r2), NaN=>0), c=:vik, clim=(-1,1))
    heatmap(replace(cor([r1 r2]), NaN=>0), c=:vik, clim=(-1,1))
    cor(sum(r1, dims=2), sum(r2, dims=2))
    plot(sum(r1, dims=2), sum(r2, dims=2), seriestype=:scatter)
end

# ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ 
# ------------------------------------------------
# Now we grab all the information for reactivation
# ------------------------------------------------
# ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑ 
coefs = OrderedDict()
iters = enumerate( zip(Z_ca1, Z_pfc))
(i, (z_ca1, z_pfc)) = first(iters)
@showprogress for (i,(z_ca1, z_pfc)) in iters
    println("Iteration:", i)
    key, v_ca1 = z_ca1
    k, v_pfc = z_pfc
    @assert(key == k)
    println("iteration", i)
    m_ca1pfc = M_PCAICA()
    m_ca1ca1 = M_PCAICA()
    m_pfcpfc = M_PCAICA()
    try
        m_ca1ca1 = reactivation.ingredients(m_ca1ca1, v_ca1, v_ca1)
        m_pfcpfc = reactivation.ingredients(m_pfcpfc, v_pfc, v_pfc)
        m_ca1pfc = reactivation.ingredients(m_ca1pfc, v_ca1, v_pfc)
    catch e
        println("Error on iteration: $s with error: $e") 
    end
    coefs[k, "ca1-pfc"] = m_ca1pfc
    coefs[k, "ca1-ca1"] = m_ca1ca1
    coefs[k, "pfc-pfc"] = m_pfcpfc
end
println("Length of coefs: ", length(coefs))


# Save the original logger
original_logger = global_logger();
# Set the global logger to suppress warnings
global_logger(ConsoleLogger(stderr, Logging.Error));
# Now we want to measure each reactivaion against all time groups,
# including those it wasn't trained on. As we do so, we will create
# a dictionary that contains those measurements taken on each.
# Before the loop ends, we will create a dataframe with the _props
# dataframe information plus the reactivation scores, and slowly
# append to it to an overall dataframe
"""
    get_df(m, z_ca1, z_pfc; kws...)
Returns a dataframe with the reactivation scores for the model `m`
on the data `z_ca1` and `z_pfc`. The keyword arguments are used
to add additional columns to the dataframe.
"""
function get_df(m, z1, z2; k_test, k_train, kws...)
    # Measure reactivation and conver to dimarray
    r = reactivation.reactscore(m, z1, z2)
    ti = Dim{:time}(props_ca1[k_test].time)
    r=DimArray(r, (ti, :component))
    # Setup dataframe from dimarray data
    df=DataFrame(r)
    df[!,:component] = convert(Vector{Int8}, df.component)
    # Append extra data from kws
    for (k,v) in kws
        df[!, Symbol(k)] .= v
    end
    # Get a dataframe of properties regarding the test and train datasets
    ptrain=DataFrame(props_ca1[k_train])
    nms = setdiff(names(ptrain),["time"])
    rename!(ptrain, nms .=> Symbol.(string.(nms, "_train")))
    ptest = DataFrame(props_ca1[k_test])
    # Register those properties to the reactivation dataframe
    for nm in setdiff(names(ptrain),["time"])
        df[!, nm] .= ptrain[1, nm]
    end
    for nm in names(ptest)
        df[!, nm] .= ptest[1, nm]
    end
    return r, df
end
ckeys = unique(map(keys(coefs) |> collect) do (prop, areas)
    prop
end)
train_iters = enumerate(Iterators.product(enumerate(ckeys), 
                                    enumerate(ckeys)))
(i,( train, test )) = collect(train_iters)[2]
Scores, DF = OrderedDict(), DataFrame()
@showprogress "Test on trainings" for (i,(train, test)) in train_iters
    i_train, k_train = train
    i_test, k_test   = test
    i, i_train, i_test = Int16.((i, i_train, i_test))
    z_ca1, z_pfc = Z_ca1[k_test], Z_pfc[k_test]
    kws = (;i, i_train, i_test)
    # Get reactivation scores
    m = coefs[k_train, "ca1-pfc"]
    r, df = get_df(m, z_ca1, z_pfc; k_train=k_train, k_test=k_test,
                   areas="ca1-pfc", kws...)
    Scores[k_train, k_test, "ca1-pfc"] = r
    append!(DF, df)
    m = coefs[k_train, "ca1-ca1"]
    r, df = get_df(m, z_ca1, z_ca1; k_train=k_train, k_test=k_test,
                   areas="ca1-ca1", kws...)
    Scores[k_train, k_test, "ca1-ca1"] = r
    append!(DF, df)
    m = coefs[k_train, "pfc-pfc"]
    r, df = get_df(m, z_pfc, z_pfc; k_train=k_train, k_test=k_test,
                   areas="pfc-pfc", kws...)
    Scores[k_train, k_test, "pfc-pfc"] = r
    append!(DF, df)
end
# Restore the original logger
global_logger(original_logger);
# Shrink memory footprint of variables
DF.i = convert(Vector{Int16}, DF.i)
DF.i_train = convert(Vector{Int16}, DF.i_train)
DF.i_test = convert(Vector{Int16},  DF.i_test)
DF.startstopWell = convert(Vector{Int8}, DF.startstopWell)
DF.startstopWell_train = convert(Vector{Int8}, DF.startstopWell_train)
DF.stopWell = convert(Vector{Int8}, DF.stopWell)
DF.startWell = convert(Vector{Int8}, DF.startWell)
DF.startWell_train = convert(Vector{Int8}, DF.startWell_train)
DF.stopWell_train = convert(Vector{Int8}, DF.stopWell_train)
DF.component = convert(Vector{Int8}, DF.component)
# Checkpoint
commit_react_vars()

# ------------------------------------------------
# Whole i_train, i_test summaries
# ------------------------------------------------
other_cols = setdiff(propertynames(DF), 
[:i, :i_train, :i_test, :component, :areas, :value, :time])
DFS = combine(groupby(DF, [:i_train, :i_test, :component, :areas]), 
        other_cols .=> [first],
        :value => mean => :mean,
        :value => std => :std,
        :value => length => :n,
        renamecols=false
)

#    _  _     _   _               _            _   _             
#  _| || |_  | | | |_   _ _ __   | |_ ___  ___| |_(_)_ __   __ _ 
# |_  ..  _| | |_| | | | | '_ \  | __/ _ \/ __| __| | '_ \ / _` |
# |_      _| |  _  | |_| | |_) | | ||  __/\__ \ |_| | | | | (_| |
#   |_||_|   |_| |_|\__, | .__/   \__\___||___/\__|_|_| |_|\__, |
#                   |___/|_|                               |___/ 
# ------------------------------------------------
DF_q = DIutils.arr.get_quantile_filtered(DF, :value, 0.001)
# Want to select train=moving and test=stationary
general_conditions = [:moving_train => x -> x .== true
                      :moving => x -> x .== false]

default_match_cols =  [[:startWell, :stopWell, :startstopWell, :ha],
[:startWell_train, :stopWell_train, :startstopWell_train, :ha_train]]
match_cols = copy(default_match_cols)
DFS.match = all(Matrix(DFS[!, match_cols[1]]) .== Matrix(DFS[!, match_cols[2]]),
                dims=2) |> vec
DFS.exclude = DFS.n .< 10


function subset_dfs(pos...)
    # Subset of DF where k_train == k_test
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

# ------------------------------------------------
# Same as above section, but split out by different startWell to the
# same stopWell given (train_startWell, train_stopWell)
#
# Step 1: Change .match to encode match excepting the startWell
# Step 2: function that subsets, hypothesis tests to the REPL, and
#         plots the histogram as above
# ------------------------------------------------
# Step 1
unmatch = :startWell
match_cols = [[x for x in default_match_cols[1] if !occursin(string(unmatch), string(x))],
              [x for x in default_match_cols[2] if !occursin(string(unmatch), string(x))]]
DFS.match = all(Matrix(DFS[!, match_cols[1]]) .== Matrix(DFS[!, match_cols[2]]),
                dims=2) |> vec;
DFS.exclude = DFS.startWell == -1 || DFS.stopWell == -1 || 
              DFS.startWell_train == -1 || DFS.stopWell_train == -1
dfs, _ = subset_dfs()
U = Dict(x=>unique(dfs[!, x]) for x in propertynames(dfs))
train_dims = [:startWell_train, :stopWell_train]
test_dims = [:startWell]
# Step 2
h = []
train_iters = Iterators.product([U[x] for x in train_dims]...)
item = first(train_iters)
for item in train_iters
    selector = [train_dims[i] => x -> x .== item[i] for i in 1:length(item)]
    dfs_match, dfs_nonmatch = subset_dfs( selector...)
    test_iters = Iterators.product([U[x] for x in test_dims]...)
    matching    = Array{Any}(undef, size(train_iters)..., size(test_iters)...)
    nonmatching = Array{Any}(undef, size(train_iters))
    startWell = first(test_iters)
    for (i,startWell) in enumerate(unique(dfs.startWell))
        dm = dfs_match[
            dfs_match.startWell .== startWell, :]
        dn = dfs_match[
            dfs_match.startWell .== startWell, :]
        if startWell_train == startWell
            push!(matching, dm)
        else
            push!(nonmatching, dn)
        end
    end
end




# ------------------------------------------------
# Register traj with DF and explore trajectory-wise reactivation scores
# ------------------------------------------------
beh.traj = replace(beh.traj, NaN=>-1)
beh.traj = convert(Vector{Int16}, beh.traj)
DIutils.filtreg.register(beh, DF, on="time", transfer=["traj"])
commit_react_vars()

