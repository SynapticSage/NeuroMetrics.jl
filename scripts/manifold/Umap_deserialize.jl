using DataStructures: OrderedDict
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
using  DimensionalData

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
using Munge.manifold

# Load data
# ----------------
@time spikes, beh, ripples, cells = Load.load("RY16", 36);
Rca1, Rpfc = (Munge.spiking.torate(@subset(spikes,:area .== ar), beh)
                for ar in ("CA1","PFC"))

# Basic params
# ----------------
filt = nothing
areas = (:ca1,:pfc)
distance = :Mahalanobis
feature_engineer = nothing
feature_engineer = :zscore
load_manis(Main; feature_engineer, filt, distance)

splits = 2
sampspersplit = 2
nsamp = Int(round(size(beh,1)/splits))
δi    = Int(round(nsamp/sampspersplit))
#nsamp = min(99_000, size(beh,1))
samps = []
for (split, samp) in Iterators.product(1:splits, 1:sampspersplit)
    start = (split-1) * nsamp + (samp-1) * δi + 1
    stop  = (split)   * nsamp + (samp)   * δi + 1
    stop  = min(stop, size(beh,1))
    push!(samps, start:stop)
end
@info "coverage" nsamp/size(beh,1)*100


keys(embedding)
key = (area = :ca1, dim = 3, s = 1, min_dist = 0.2, n_neighbors = 200)
em = embedding[key]
em = DimArray(em, (Dim{:time}(beh.time[1:size(em,1)]), Dim{:comp}(1:size(em,2))))

