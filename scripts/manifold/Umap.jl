using DataStructures: OrderedDict
using DrWatson
using  GoalFetchAnalysis
import Utils.namedtup: ntopt_string
using  DataFramesMeta
using  Serialization
using Plots
using ProgressMeter
using PyCall
using Distributed
using ProgressMeter
using ThreadSafeDicts
using DataFramesMeta
using Distances
using StatsBase

# Disstributed computing
#addprocs([("mingxin",1)])
#addprocs(2)
#pids = workers()
@everywhere using DrWatson
@everywhere quickactivate(expanduser("~/Projects/goal-code"))
@everywhere using PyCall
@everywhere using  UMAP
@everywhere using DataFramesMeta
@everywhere using GoalFetchAnalysis
@everywhere import Munge

# Load data
# ----------------
@time spikes, beh, ripples, cells = Load.load("RY16", 36);

# Basic params
# ----------------
filt = nothing
areas = (:ca1,:pfc)
distance = :Mahalanobis
feature_engineer = :zscore
use_cuda = true
if use_cuda
    ENV["PYTHON"]="/home/ryoung/miniconda3/envs/rapids-22.08/bin/python"
    using PyCall
    cuml = pyimport("cuml")
end

# Filter
if filt !== nothing
    filtstr = "filt=$filt"
    filters = Filt.get_filters()[filt]
    Utils.filtreg.filterAndRegister(beh, spikes; filters)
else
    filtstr = "filt=nothing"
end
festr   = feature_engineer === nothing ? "feature=nothing" : "feature=$feature_engineer"
diststr = distance === nothing ? "distance=euclidean" : lowercase("distance=$distance")
savefile = datadir("manifold","ca1pfc_manifolds_$(filtstr)_$(diststr)_$(festr).serial")
@info "run info" filtstr festr diststr savefile

# Get sample runs
# ----------------
splits, sampspersplit = 3, 2
nsamp = Int(round(size(beh,1)/splits))
δi    = Int(round(nsamp/sampspersplit))
samps = []
for (split, samp) in Iterators.product(1:splits, 1:sampspersplit)
    start = (split-1) * nsamp + (samp-1) * δi + 1
    stop  = (split)   * nsamp + (samp-1) * δi + 1
    stop  = min(stop, size(beh,1))
    push!(samps, start:stop)
end
@info "coverage" nsamp/size(beh,1)*100

# Get rate matrices
# ----------------
R  = Dict()
R[:ca1], R[:pfc] = (Munge.spiking.torate(@subset(spikes,:area .== ar), beh)
                    for ar in ("CA1","PFC"))
if feature_engineer == :zscore
    for area in areas
        R[area ] = hcat([zscore(c) for c  in eachcol(R[area])]...)
    end
end
if distance == :Mahalanobis
    Q = Dict()
    dist_func = Dict()
    for area in areas
        Q[area] = R[area]' * R[area]
        dist_func[area] = getproperty(Distances, distance)(Matrix(Q[area]))
    end
else
    dist_func = Dict(area=>getproperty(Distances, distance)
                     for area in areas)
end

# Get embeddings
# ----------------
embedding    = Dict()
dimset       = (2,   3)
min_dists    = (0.3, 0.3)
n_neighborss = (150, 200)

(min_dist, n_neighbors) = first(zip(min_dists, n_neighborss))
(area,dim,s) = first(Iterators.product(areas, dimset, (1,)))

@showprogress "params" for (min_dist, n_neighbors) in zip(min_dists, n_neighborss)
    @showprogress "datasets" for (area,dim,s) in Iterators.product(areas, dimset, (1,))
        key = (;area,dim,s,min_dist,n_neighbors)
        if key ∈ keys(embedding) && !(embedding[key] isa Future)
            @info "skipping" key
            continue
        else
            @info key
        end
        input = Matrix(R[area]'[:,samps[s]])
        em = if use_cuda
            nothing
        else
            umap(input, dim; min_dist, n_neighbors, 
                        metric=dist_func[area])
        end
        embedding[key] = em'
    end
end

using SoftGlobalScope
#embedding = Dict(k=>(if isready(v); fetch(v); else; v; end) for (k,v) in embedding); 
embedding = Dict(k=>(try; fetch(v); catch; v; end) for (k,v) in embedding);
#display(embedding); 


# ----------------------------------------------------
# Get clean indices that are within quantile tolerance
# ----------------------------------------------------

# Store them for later
@info "save info" filtstr festr diststr savefile
serialize(savefile, (;embedding, inds, samps, animal="RY16", day=36))

using Serialization
embedding, inds, animal, day = deserialize(savefile)
