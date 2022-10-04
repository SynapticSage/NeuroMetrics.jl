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
filt = :all
areas = (:ca1,:pfc)
distance = :Mahalanobis
feature_engineer = :zscore

# Filter
if filt !== nothing
    filtstr = "_$filt"
    filters = Filt.get_filters()[filt]
    Utils.filtreg.filterAndRegister(beh, spikes; filters)
else
    filtstr = ""
end
festr   = feature_engineer === nothing ? "" : "feature=$feature_engineer"
diststr = distance === nothing ? "" : "distance=$feature_engineer"
@info filtstr

# Get sample runs
# ----------------
splits = 3
sampspersplit = 3
nsamp = Int(round(size(beh,1)/splits))
δi    = Int(round(nsamp/sampspersplit))
#nsamp = min(100_000, size(beh,1))
samps = []
for (split, samp) in Iterators.product(1:splits, 1:sampspersplit)
    start = (split-1) * nsamp + (samp-1) * δi + 1
    stop  = (split)   * nsamp + (samp)   * δi + 1
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

dimset       = (2,3)
min_dists    = (0.5,0.1, 0.2,0.5)
n_neighborss = (30,15, 200, 200)

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
        embedding[key] =
        umap(Matrix(R[area]'[:,samps[s]]), dim; min_dist, n_neighbors, 
                    metric=dist_func[area])'
    end
end

using SoftGlobalScope
#embedding = Dict(k=>(if isready(v); fetch(v); else; v; end) for (k,v) in embedding); 
embedding = Dict(k=>(try; fetch(v); catch; v; end) for (k,v) in embedding);
#display(embedding); 


# Get clean indices that are within quantile tolerance
# ----------------
inds = Dict()
for key in keys(inds)
    inds[key] = Utils.clean.inds_quantile_filter_dims(
                        embedding[key], [0.02, 0.96])
end

# Store them for later
serialize(datadir("manifold","ca1pfc_manifolds$(filtstr)$(diststr)$(festr).serial"),
          (;embedding, inds, animal="RY16", day=36))

using Serialization
embedding, inds, animal, day = deserialize(datadir("manifold","ca1pfc_manifolds$(filtstr).serial"))

