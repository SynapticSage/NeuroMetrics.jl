if defined(Main, :DrWatson)
    include(scriptsdir("reactivate","imports.jl"))
else
    include("imports.jl") # include(scriptsdir("reactivate","imports.jl"))
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
println("Available trajectories:", unique(beh[:, :traj]))

# Register that to Rdf
sort!(Rdf, [:time, :unit])
register(beh, Rdf,   on="time", transfer=["traj",string.(groupings)...])
Rdf = @subset(Rdf, :stopWell .!= -1, :startWell .!= 1)
sort!(Rdf, [:unit, :time])
register(cells, Rdf, on="unit", transfer=["area"])
sort!(Rdf, [:time, :unit])
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
# INFO: up to this point, all startWell elements exist!
println("Available trajectories, Rdf:", unique(Rdf[:, :traj]))
println("Available trajectories, ca1:", unique(combine(ca1,identity)[:, :traj]))

# IF startWell=1 is actually, being rejected, I need to actually plot the conditions
# that lead to its rejection
@subset(Rdf, :startWell .== 1) |> Voyager()
@subset(beh, :startWell .== 1) |> Voyager()

GC.gc()

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
    # This line works as follows: 
    # 1. Take the first row of each cell group
    # 2. Filter out the groups that have less than num_cells samples
    # 3. Filter out the groups that have less than 3 trajectories
    # 4. Map the groups to a new dataframe with the groupings as columns
    #    and the size of the group and number of unique trajectories in the
    #    group as new columns
    accepted = 
        Base.map(
        Base.filter(map(x->groupby(x,:unit)|>first ,collect(ca1))) do x #1
        size(x,1) > size(@subset(cells,:area.=="CA1"),1) && length(unique(x.traj)) ≥ 3 #2,3
    end) do x #4
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
# INFO: `accepted` :: there are no accepted startWell=1 combinations for RY16,36


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
println("Available trajectories, ca1:", 
    unique(combine(ca1,identity)[:, :traj]))
println("Number of trajectories, ca1:", 
    length(unique(combine(ca1,identity)[:, :traj])))

# ----------------------------------
# Now convert these back to DimArrays
# ----------------------------------
ca1 = map(zip(keys(ca1),ca1) |> collect) do (k,x)
    NamedTuple(k) => unstack(x, [:time,:traj,groupings...], 
        :unit, :rate; combine=mean)
end
pfc = map(zip(keys(pfc),pfc) |> collect) do (k,x)
    NamedTuple(k) => unstack(x, [:time,:traj,groupings...], 
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
props = [groupings..., :time, :traj]
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
    zca1, zpfc = Z_ca1 |> values |> first, Z_pfc |> values |> first
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
    plot(sum(r1, dims=2), sum(r2, dims=2), seriestype=:scatter,
        xlabel="CA1-CA1", ylabel="CA1-PFC", legend=false
    )
    # Normalize each component to be length 1
    # (this might better be served by projection to search for
    # entangled and disentangled components)
    P, Curling, Div, DC = [], [], [], []
    for k in 1:min(size(r1,2),size(r2,2))
        rr1 = r1[:,1:k]
        rr2 = r2[:,1:k]
        rr1 = rr1 ./ sqrt.(sum(rr1.^2, dims=1)) ./ size(rr1,2)
        rr2 = rr2 ./ sqrt.(sum(rr2.^2, dims=1)) ./ size(rr2,2)
        rr1 = nansum(rr1, dims=2)
        rr2 = nansum(rr2, dims=2)
        p=plot(rr1, rr2, seriestype=:scatter,
            xlabel="CA1-CA1", ylabel="CA1-PFC", legend=false
        )
        curling = [rr1 rr2] * [-1; 1] # there may be a better way to get the
                                      # curling component
        h=histogram(curling, bins=100, xlabel="Curling", ylabel="Count",
            edgealpha=0, linewidth=0, fillalpha=0.5, legend=false
        )
        vline!([0], c=:black, linewidth=2, linestyle=:dash, legend=false,
        )
        diving = [rr1 rr2] * [1; 1]
        h2=histogram(diving, bins=100, xlabel="Diving", ylabel="Count",
            edgealpha=0, linewidth=0, fillalpha=0.5, legend=false
        )
        h3=plot(
            # histogram(abs.(diving)./abs.(curling), bins=100, xlabel="Diving/Curling", ylabel="Count", edgealpha=0, linewidth=0, fillalpha=0.5, legend=false), 
            plot(abs.(diving), abs.(curling), seriestype=:scatter,
                    xlabel="Diving", ylabel="Curling", legend=false
                ))
        push!(P, p)
        push!(Curling, h)
        push!(Div, h2)
        push!(DC, h3)
    end
    plot(P..., marker=:circle, markersize=1, fontsize=3, tickfontsize=3,
        labelfontsize=3, legendfontsize=3)
    plot(Curling..., marker=:circle, markersize=1, fontsize=3, tickfontsize=3,
        labelfontsize=3, legendfontsize=3)
    plot(Div..., marker=:circle, markersize=1, fontsize=3, tickfontsize=3,
        labelfontsize=3, legendfontsize=3)
    plot(DC..., marker=:circle, markersize=1, fontsize=3, tickfontsize=3, fillalpha=0.2, linewidth=0, get_markerstrokewidth=0, markerstrokealpha=0, labelfontsize=3, legendfontsize=3)
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


# ------------------------------------------------
#             MUNGE to a dataframe
# ------------------------------------------------
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
function get_df(m, z1, z2; k_test, k_tmpl, kws...)
    # Measure reactivation and conver to dimarray
    r = reactivation.reactscore(m, z1, z2)
    ti = Dim{:time}(props_ca1[k_test].time)
    r  = DimArray(r, (ti, :component))
    # Setup dataframe from dimarray data
    df=DataFrame(r)
    df[!,:component] = convert(Vector{Int8}, df.component)
    # Append extra data from kws
    for (k,v) in kws
        df[!, Symbol(k)] .= v
    end
    # Get a dataframe of properties regarding the test and train datasets
    ptrain = DataFrame(props_ca1[k_tmpl])
    nms    = setdiff(names(ptrain),["time"])
    rename!(ptrain, nms .=> Symbol.(string.(nms, "_tmpl")))
    ptest = DataFrame(props_ca1[k_test])
    # Register those properties to the reactivation dataframe
    dropmissing!(ptrain)
    for nm in setdiff(names(ptrain),["time"])
        df[!, nm] .= ptrain[1, nm]
    end
    dropmissing!(ptest)
    for nm in setdiff(names(ptest),["time"])
        df[!, nm] .= ptest[1, nm]
    end
    sort!(df, [:time, :component])
    DIutils.filtreg.register(ptest, df, on="time", transfer=["traj"]) 
    return r, df
end
# ISSUE: Already missing my startWell=1
ckeys = unique(map(keys(coefs) |> collect) do (prop, areas)
    prop
end)
train_iters = enumerate(Iterators.product(enumerate(ckeys), 
              enumerate(ckeys)))
(i,( train, test )) = collect(train_iters)[2]
Scores, DF = OrderedDict(), DataFrame()
@showprogress "Test on trainings" for (i,(train, test)) in train_iters
    i_tmpl, k_tmpl = train
    i_test, k_test   = test
    i, i_tmpl, i_test = Int16.((i, i_tmpl, i_test))
    z_ca1, z_pfc = Z_ca1[k_test], Z_pfc[k_test]
    kws = (;i, i_tmpl, i_test)
    # Get reactivation scores
    m = coefs[k_tmpl, "ca1-pfc"]
    r, df = get_df(m, z_ca1, z_pfc; k_tmpl=k_tmpl, k_test=k_test,
                   areas="ca1-pfc", kws...)
    Scores[k_tmpl, k_test, "ca1-pfc"] = r
    append!(DF, df)
    m = coefs[k_tmpl, "ca1-ca1"]
    r, df = get_df(m, z_ca1, z_ca1; k_tmpl=k_tmpl, k_test=k_test,
                   areas="ca1-ca1", kws...)
    Scores[k_tmpl, k_test, "ca1-ca1"] = r
    append!(DF, df)
    m = coefs[k_tmpl, "pfc-pfc"]
    r, df = get_df(m, z_pfc, z_pfc; k_tmpl=k_tmpl, k_test=k_test,
                   areas="pfc-pfc", kws...)
    Scores[k_tmpl, k_test, "pfc-pfc"] = r
    append!(DF, df)
    if i == 1
        DF.i                  = convert(Vector{Int16}, DF.i)
        DF.i_tmpl             = convert(Vector{Int16}, DF.i_tmpl)
        DF.i_test             = convert(Vector{Int16}, DF.i_test)
        DF.startstopWell      = convert(Vector{Int8}, DF.startstopWell)
        DF.startstopWell_tmpl = convert(Vector{Int8}, DF.startstopWell_tmpl)
        DF.stopWell           = convert(Vector{Int8}, DF.stopWell)
        DF.startWell          = convert(Vector{Int8}, DF.startWell)
        DF.startWell_tmpl     = convert(Vector{Int8}, DF.startWell_tmpl)
        DF.stopWell_tmpl      = convert(Vector{Int8}, DF.stopWell_tmpl)
        DF.component          = convert(Vector{Int8}, DF.component)
    end
end
# Restore the original logger
global_logger(original_logger);
# Shrink memory footprint of variables
sort!(DF, [:time, :i])
# Rename any vars with _train to _tmpl
nms = names(DF)[occursin.("_train", names(DF))]
rename!(DF, nms .=> replace.(nms, "_train" => "_tmpl"))
println("Number of trajectories: $(length(unique(DF.traj)))")

# ------------------------------------------------
# Register traj with DF and explore trajectory-wise reactivation scores
# ------------------------------------------------
# beh.traj = replace(beh.traj, NaN=>-1)
# beh.traj = convert(Vector{Int16}, beh.traj)
# DIutils.filtreg.register(beh, DF, on="time", transfer=["traj"])

# ------------------------------------------------
# Whole i_tmpl, i_test summaries
# ------------------------------------------------
other_cols = setdiff(propertynames(DF), 
[:i, :i_tmpl, :i_test, :component, :areas, :value, :time])
DFS = combine(groupby(DF, [:i_tmpl, :i_test, :component, :areas]), 
        other_cols .=> [first],
        :value => mean => :mean,
        :value => std => :std,
        :value => length => :n,
        renamecols=false
)

# Checkpoint
commit_react_vars()

