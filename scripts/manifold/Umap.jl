quickactivate(expanduser("~/Projects/goal-code"))
using DataStructures: OrderedDict
using DrWatson
using Serialization
using Plots
using ProgressMeter
using PyCall
using Distributed
using ProgressMeter
using ThreadSafeDicts
using DataFramesMeta
using Distances
using StatsBase
use_cuda = true
if use_cuda
    #ENV["PYTHON"]="/home/ryoung/miniconda3/envs/rapids-22.08/bin/python"
    using PyCall
    cuml = pyimport("cuml")
    cuUMAP = cuml.manifold.umap.UMAP
end

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
import Utils.namedtup: ntopt_string

# Load data
# ----------------
@time spikes, beh, ripples, cells = Load.load("RY16", 36)

# Basic params
# ----------------
filt             = nothing
areas            = (:ca1,:pfc)
#distance        = :Mahalanobis
distance         = :CityBlock
feature_engineer = :zscore

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
embedding,scores    = Dict(),Dict()
dimset       = (2,   3)
min_dists    = (0.05,0.15,0.3)
n_neighborss = (5,50,150,400)

(min_dist, n_neighbors) = first(zip(min_dists, n_neighborss))
(area,dim,s) = first(Iterators.product(areas, dimset, (1,)))

@showprogress "params" for (min_dist, n_neighbors) in Iterators.product(min_dists, n_neighborss)
    @showprogress "datasets" for (area,dim,s) in Iterators.product(areas, dimset, 
                                                                   1:length(samps))
        key = (;area,dim,s,min_dist,n_neighbors)
        if key ∈ keys(embedding) && !(embedding[key] isa Future)
            @info "skipping" key
            continue
        else
            @info key
        end
        input = Matrix(R[area]'[:,samps[s]])
        @time em, score = if use_cuda # 1000x faster
            metric = lowercase(String(distance))
            fitter=cuUMAP(n_neighbors=n_neighbors, min_dist=min_dist, 
                          n_components=3, metric=metric, 
                          target_metric="euclidean")
            trained_umap = fitter.fit(input')
            em = trained_umap.transform(input')
            score = cuml.metrics.trustworthiness(input', em, 
                                                 n_neighbors=n_neighbors)
            em, score
        else # julia is suprisingly slow here
            em = umap(input, dim; min_dist, n_neighbors, 
                        metric=dist_func[area])
            skmani = pyimport("sklearn.manifold")
            @time score = skmani.trustworthiness(input', em, 
                                           n_neighbors=n_neighbors, 
                                           metric=lowercase(string(distance)))
            em, score

        end
        embedding[key] = em'
        scores[key] = score
    end
end

using SoftGlobalScope
embedding = Dict(k=>(try; fetch(v); catch; v; end) for (k,v) in embedding);

#embedding = Dict(k=>(if isready(v); fetch(v); else; v; end) for (k,v) in embedding); 
#display(embedding); 


# ----------------------------------------------------
# Get clean indices that are within quantile tolerance
# ----------------------------------------------------

# Store them for later
@info "save info" filtstr festr diststr savefile
save_manis(;embedding, filt, feature_engineer, distance, samps, use_cuda)

using Serialization
embedding, inds, animal, day = deserialize(savefile)
