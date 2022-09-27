using DataStructures: OrderedDict

using DrWatson
using  GoalFetchAnalysis
import Utils.namedtup: ntopt_string
using  DataFramesMeta
using  UMAP
using  Serialization
using Plots
using ProgressMeter
using PyCall
using ThreadSafeDicts

# Load data
# ----------------
@time spikes, beh, ripples, cells = Load.load("RY16", 36);

# Basic params
# ----------------
filt = :all
if filt !== nothing
    filtstr = "_$filt"
    filters = Filt.get_filters()[filt]
    Utils.filtreg.filterAndRegister(beh, spikes; filters)
else
    filtstr = ""
end
@info filtstr
areas = (:ca1,:pfc)

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

# Get embeddings
# ----------------
embedding    = Dict()
dimset       = (2,3,11)
min_dists    = (0.5,0.1, 0.2)
n_neighborss = (30,15, 200)
@showprogress "params" for (min_dist, n_neighbors) in zip(min_dists, n_neighborss)
    @showprogress "datasets" for (area,dim,s) in Iterators.product(areas,dimset, (1,))
        embedding[(area,dimset,s)]  = umap((R[area]')[:,samps[s]], dim; min_dist, n_neighbors)'
    end
end

# Get clean indices that are within quantile tolerance
# ----------------
inds = Dict()
for (area,dim) in Iterators.product(areas,dimset)
    inds[(area,dim)] = Utils.clean.inds_quantile_filter_dims(
                        embedding[area], [0.02, 0.96])
end

# Store them for later
serialize(datadir("manifold","ca1pfc_manifolds$(filtstr).serial"),
          (;embedding, inds, animal="RY16", day=36))


