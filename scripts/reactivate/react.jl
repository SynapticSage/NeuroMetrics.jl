include("imports.jl")
if !isdefined(Main, :ca1) && !isdefined(Main, :pfc) && !isdefined(Main, :opt)
    include("run.jl")
end

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
    plot(DC..., marker=:circle, markersize=1, fontsize=3, tickfontsize=3, fillalpha=0.2, linewidth=0, 
        markerstrokewidth=0, markerstrokealpha=0, labelfontsize=3, legendfontsize=3)
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
    df[!,:component] = convert(Vector{UInt8}, df.component)
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
    # These steps are essential BUT they are costly ... on my machine
    # adds 15 minutes to a 20 minute loop.
        sort!(df, [:time, :component])
        DIutils.filtreg.register(ptest, df, on="time", transfer=["traj"]) 
    return r, df
end
function convert_df!(df)
    df.i                  = convert(Vector{Int16}, df.i)
    df.i_tmpl             = convert(Vector{Int16}, df.i_tmpl)
    df.i_test             = convert(Vector{Int16}, df.i_test)
    df.startstopWell      = convert(Vector{Int8},  df.startstopWell)
    df.startstopWell_tmpl = convert(Vector{Int8},  df.startstopWell_tmpl)
    df.stopWell           = convert(Vector{Int8},  df.stopWell)
    df.startWell          = convert(Vector{Int8},  df.startWell)
    df.startWell_tmpl     = convert(Vector{Int8},  df.startWell_tmpl)
    df.stopWell_tmpl      = convert(Vector{Int8},  df.stopWell_tmpl)
    df.component          = convert(Vector{UInt8},  df.component)
    df
end
# ISSUE: Already missing my startWell=1
ckeys = unique(map(keys(coefs) |> collect) do (prop, areas)
    prop
end)
train_iters = enumerate(Iterators.product(enumerate(ckeys), 
              enumerate(ckeys))) |> collect
(i,( train, test )) = collect(train_iters)[2]
prog = Progress(length(train_iters); desc="Reactivation Scores")
Scores, DF = OrderedDict(), Matrix{Union{Missing,DataFrame}}(missing, length(train_iters), 3)
println("Garbage colleting:", GC.gc())
for (i,(train, test)) in train_iters
    i_tmpl, k_tmpl   = train
    i_test, k_test   = test
    i, i_tmpl, i_test = Int16.((i, i_tmpl, i_test))
    z_ca1, z_pfc = Z_ca1[k_test], Z_pfc[k_test]
    kws = (;i, i_tmpl, i_test)
    # Get reactivation scores
    m = coefs[k_tmpl, "ca1-pfc"]
    r, df = get_df(m, z_ca1, z_pfc; k_tmpl=k_tmpl, k_test=k_test,
                   areas="ca1-pfc", kws...)
    Scores[k_tmpl, k_test, "ca1-pfc"] = r
    DF[i, 1] = convert_df!(df)
    m = coefs[k_tmpl, "ca1-ca1"]
    r, df = get_df(m, z_ca1, z_ca1; k_tmpl=k_tmpl, k_test=k_test,
                   areas="ca1-ca1", kws...)
    Scores[k_tmpl, k_test, "ca1-ca1"] = r
    DF[i, 2] = convert_df!(df)
    m = coefs[k_tmpl, "pfc-pfc"]
    r, df = get_df(m, z_pfc, z_pfc; k_tmpl=k_tmpl, k_test=k_test,
                   areas="pfc-pfc", kws...)
    Scores[k_tmpl, k_test, "pfc-pfc"] = r
    DF[i, 3] = convert_df!(df)
    next!(prog)
end
# Restore the original logger
global_logger(original_logger);

# Need to now combine the dataframes
DFnew = DataFrame()
DF = DF[:]
while !isempty(DF)
    append!(DFnew, popfirst!(DF))
end
if !isempty(DFnew)
    DF = DFnew
end
GC.gc()


sort!(DF, [:time, :i])
# Rename any vars with _train to _tmpl
nms = names(DF)[occursin.("_train", names(DF))]
rename!(DF, nms .=> replace.(nms, "_train" => "_tmpl"))
println("Number of trajectories: $(length(unique(DF.traj)))")

# Ensure that .component is a UInt8 and that each template has a normal range
# of trajectories
tmp = combine(groupby(DF, [:startWell_tmpl,  :stopWell_tmpl]),
    :component => (x->length(unique(x))) =>  :count,
    :component => (x->maximum(unique(x))) => :max,
    :component => (x->minimum(unique(x))) => :min
)
tmp.startstop = string.(tmp.startWell_tmpl, "-", tmp.stopWell_tmpl)
@df tmp bar(:startstop, :count, title="Number of components per template")
@df tmp scatter(:startstop, :max, title="Number of components per template")
@df tmp scatter!(:startstop, :min, title="Number of components per template")
# Draw a line between the points
p = plot!()
plot!([p.series_list[1][:x] p.series_list[1][:x]]', [tmp.min tmp.max]', 
    title="Number of components per template", label="", color=:black)

# ------------------------------------------------
# Register traj with DF and explore trajectory-wise reactivation scores
# ------------------------------------------------
# beh.traj = replace(beh.traj, NaN=>-1)
# beh.traj = convert(Vector{Int16}, beh.traj)
# DIutils.filtreg.register(beh, DF, on="time", transfer=["traj"])

# ------------------------------------------------
# Create a time-averaged summary DF
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

# ------------------------------------------------
# Create a component averaged DF
# ------------------------------------------------
# function Statistics.mean(x::AbstractArray{Union{Missing, String}})
#     x = x |> skipmissing |> collect
#     @assert all(x .== x[1])
#     first(x)
# end
# function Statistics.median(x::AbstractArray{Union{Missing, String}})
#     x = x |> skipmissing |> collect
#     @assert all(x .== x[1])
#     first(x)
# end
nms = setdiff(names(DF), ["component","value"])
tmp = groupby(subset(DF, view=true), nms)
DF1 = combine(
    tmp,
    :value => mean   => :mean,
    :value => std    => :std,
    :value => length => :n,
    :value => median => :median,
    renamecols=false
)
println("Size fraction: $(round(size(DF1,1)./size(DF,1), digits=4))")

# Checkpoint
@time commit_react_vars()
