#!/bin/sh
#=
export LD_LIBRARY_PATH=/home/ryoung/miniconda3/envs/conda_jl/lib/
exec julia -J "/home/ryoung/Code/projects/goal-code/GFA-dependencies-sysimage.so" --project="/home/ryoung/Projects/goal-code/" "$0" -- $@
=#
using GoalFetchAnalysis, Munge.manifold, DrWatson, Revise
opt = Munge.manifold.parse()

# Sets to explore
# ---------------
#metrics      = unique((:CityBlock, :Euclidean,:Correlation,:Cosine))
#dimset       = (2,   3)
#min_dists    = (0.05,0.15,0.3)
#n_neighborss = (5,50,150,400)
min_dists, n_neighborss, metrics, dimset, features = [0.3], [5,50], [:Euclidean], 
                                                     [2,3], [:zscore]


# FINISH loading libraries
begin
    using Plots: StatsBase
    using Infiltrator, Serialization, Plots, ProgressMeter, PyCall, Distributed,
          ProgressMeter, ThreadSafeDicts, DataFramesMeta, Distances, StatsBase,
          SoftGlobalScope, Infiltrator, DimensionalData,  DataFramesMeta
    using Plot.manifold
    using DataStructures: OrderedDict
    import DIutils.namedtup: ntopt_string
    use_cuda = true
    if use_cuda
        using PyCall
        cuml = pyimport("cuml")
        gc = pyimport("gc")
        cuUMAP = cuml.manifold.umap.UMAP
    end
end

for (key, value) in opt
    if !(value isa Symbol)
        @eval Main global $(Symbol(key)) = $value
    else
        value = String(value)
        @eval Main global $(Symbol(key)) = Symbol($value)
    end
end
#global areas            = (:ca1,:pfc)

# Load data
# ----------------


println("Loading")
@time global spikes, beh, ripples, cells = Load.load(opt["animal"], opt["day"])
cells, spikes = DIutils.filtreg.register(cells, spikes; on="unit", transfer=["celltype"])
beh.index = 1:size(beh,1)

# ----------------
# Basic params
# ----------------
global spikesTrain = spikesTest = spikes
global behTrain = behTest = beh

# Filter
println("Filtration?")
filt, distance, feature_engineer = opt["filt"], opt["distance"], opt["feature_engineer"]
if filt !== nothing
    global filtstr = "filt=$filt"
    filters = Filt.get_filters()[filt]
    global behTrain, spikesTrain =
                        DIutils.filtreg.filterAndRegister(copy(behTrain), copy(spikesTrain); filters,
                                                       filter_skipmissingcols=true)
else
    global filtstr = "filt=nothing"
end
global festr   = opt["feature_engineer"] === nothing ? "feature=nothing" : "feature=$feature_engineer"
global diststr = opt["distance"] === nothing ? "distance=euclidean" : lowercase("distance=$distance")
@info "run info" filtstr festr diststr 

# Firing Rates
# ------------
println("Firing rate matrices")
function get_R(beh, spikes)
    
    R = Dict(
             Symbol(lowercase(ar)) =>
             Munge.spiking.torate(@subset(spikes,:area .== ar), beh)
                    for ar in ("CA1","PFC")
    )
    R = merge(R, 
              Dict(
                  Symbol(lowercase(ar) * "_" * String(ct)) => 
                  (sub=@subset(spikes, :area .== ar, :celltype .== ct);
                   if !isempty(sub)
                       Munge.spiking.torate(sub, beh, gaussian=4)
                   else
                       []
                   end
                  ) for ar in ("CA1","PFC"), ct in (:pyr,:int)
              )
    )
    R = Dict(k=>v for (k,v) in R if v != [])
    zscoredimarray(x) = DimArray(hcat(zscore.(eachcol(x))...), x.dims)
    merge(R,Dict(Symbol("z"*string(k))=>zscoredimarray(v) for (k,v) in R if v != []))
end
begin
    R      = get_R(beh, spikes)
    # Rtrain = get_R(behTrain, spikesTrain)
    Rtrain_inds = Dict(k=>DIutils.searchsortednearest.([beh.time], behTrain.time) for (k,_) in R)
    Rtrain = Dict(
                  k=>v[inds,:] for ((k,v),(k,inds)) in zip(R,Rtrain_inds)
                 )
end

# Get sample runs
# ----------------
println("Generate partitions")
splits, sps = opt["splits"], opt["sps"]
nsamp = Int(round(size(beh,1)/splits));
δi    = Int(round(nsamp/sps))
global inds_of_t = []
for (split, samp) in Iterators.product(1:splits, 1:sps)
    start = (split-1) * nsamp + (samp-1) * δi + 1
    stop  = (split)   * nsamp + (samp-1) * δi + 1
    stop  = min(stop, size(beh,1))
    push!(inds_of_t, start:stop)
end

println("Describ partitions")
@info "coverage" nsamp/size(beh,1)*100
global N = opt["splits"] * opt["sps"]

global tag 
tag = "$(opt["animal"])$(opt["day"]).$(N)seg"
println(tag)

# Find the indices in the training sets that match each partition
Rtrain_matching = Dict(
           (k, inds) => begin
               q(i) = i >= first(inds) && i <= last(inds)
               findall(q.(v))
               
           end
           for (k,v) in Rtrain_inds, inds in inds_of_t
)

# ----------------
# Get prior embeddings
# ----------------
global embedding, scores = if isfile(path_manis(;filt,feature_engineer,tag))
    @info "loading prev data"
    data=load_manis(Main;filt,feature_engineer,tag);
    @assert(all(data[:inds_of_t] .== inds_of_t))
    embedding, scores = data.embedding, data.scores
else
    embedding, scores = Dict(), Dict()
end
@info "length of dicts" length(embedding) length(scores)

# ----------------
# Setup loop variables and progress measures
# ----------------
params   = collect(Iterators.product(metrics,min_dists, n_neighborss,features))
datasets = collect(Iterators.product(keys(R), dimset, 1:length(inds_of_t)))
prog = Progress(prod(length.((params,datasets))); desc="creating embeddings")
global steps, total = 0, (length(params) * length(datasets))

trained_umap = em = sc = nothing
exception_triggered = false

(metric,min_dist, n_neighbors, feature) = first(params)
(dataset,dim,s) = first(datasets)

gc.collect()
GC.gc(false)
gc.collect()
GC.gc(false)

dim         = 3
n_neighbors = 400
dataset     = :zpfc
min_dist = 0.5
key = (;dataset,dim,s,min_dist,n_neighbors,metric,feature)

# Pre - process : Do we obtain a distance function?
@debug "Getting distance metric"
metric_str = if distance == :Mahalanobis
    @info "transforming Mahalanobis"
    Q = Rtrain[dataset]' * Rtrain[dataset];
    dist_func = getproperty(Distances, distance)(Matrix(Q[dataset]));
elseif !use_cuda
    @info "transforming $distance"
    dist_func = getproperty(Distances, distance)
else
    lowercase(string(metric))
end

train = Matrix(Rtrain[dataset]'); # train on everything with filters
#train = Matrix(Rtrain[area]'[:, # train on just this set of indices (dynamic manifold)
                             # Rtrain_matching[area, inds_of_t[s]]
                            # ]); 
# input = Matrix(R[dataset]'[:, # test on this set minus filters
#                         inds_of_t[s]
#                        ]);
input = Matrix(R[dataset]'); # test on everything minus filters

@debug "Processing"
@debug "fitter"

fitter=cuUMAP(n_neighbors=n_neighbors, min_dist=min_dist, n_components=dim, metric=metric_str, local_connectivity=1, target_metric="euclidean", n_epochs=1_000);

# local_connectivity: int (optional, default 1)
#     The local connectivity required -- i.e. the number of nearest
#     neighbors that should be assumed to be connected at a local level.
#     The higher this value the more connected the manifold becomes
#     locally. In practice this should be not more than the local intrinsic
#     dimension of the manifold.
# n_epochs larger = more accurate, 500 for small, 200 for large
@debug "fit"
T = train
trained_umap = fitter.fit(T');
@debug "transform"
# randsamp(x,n) = x[:,Random.randperm(size(x,2))[1:n]]
# I = randsamp(input, 50_000)
I = input
em = trained_umap.transform(I');
@debug "trust"
#sc = cuml.metrics.trustworthiness(I', em, n_neighbors=n_neighbors)
# t=@task Plot.stereoscopicgif( eachcol(E)... ;deltaangle=5, xlims=(-w,w), ylim=(-w,w), zlim=(-w,w), alphaasync =0.1*k)



#    _  _     ____  _          _   
#  _| || |_  |  _ \| |    ___ | |_ 
# |_  ..  _| | |_) | |   / _ \| __|
# |_      _| |  __/| |__| (_) | |_ 
#   |_||_|   |_|   |_____\___/ \__|
#                                  
using Plot.manifold
maniplot(em; by=beh.stopWell, llim=30)

using Blink, Interact
if isdefined(Main, :w); close(w); end
M = Dict()
@showprogress for θ=0:5:360, ρ=0:5:360, origin=0:10:100
    🔑 = (θ, ρ, origin)
    M[🔑] = 🔑 ∉ keys(M) ? Plot.manifold.maniplot(em; by=beh.cuemem, size=(1200,1200),
                                llim=30, camera=(θ,ρ), origin) : 
                            M[🔑]
end
ui= @manipulate for θ=0:5:360, ρ=0:5:360, origin=0:10:100
    🔑 = (θ, ρ, origin)
    M[🔑]
end
w=Window()
body!(w, ui)



B = groupby(beh, [:epoch, :period])
P = []
for (i,b) in enumerate(B)
    @info i
    title="$(b.startWell[1]), $(b.stopWell[1])"
    push!(P,
          maniplot(em[b.index,:];
                   by=b.time, pointalpha=0.5, 
                   linealpha=0.5, 
                   random=false,
                   byaxis=true,
                   cmap=:vik, llim=30, title)
    )
end
plot(P[1:10]...)
        
