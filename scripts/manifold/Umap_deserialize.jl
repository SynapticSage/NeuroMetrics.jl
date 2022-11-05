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
using SoftGlobalScope
use_cuda = true
#if use_cuda
#    #ENV["PYTHON"]="/home/ryoung/miniconda3/envs/rapids-22.08/bin/python"
#    using PyCall
#    cuml = pyimport("cuml")
#    cuUMAP = cuml.manifold.umap.UMAP
#end
using DimensionalData

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
datasets = ( ("RY16", 36),)
filt = nothing
areas = (:ca1,:pfc)
distance = :many
feature_engineer = nothing
feature_engineer = :many
N = 100
embedding_overall = Dict()
@softscope for (animal,day) in datasets
    
    #@time spikes, beh, ripples, cells = Load.load(animal,day)
    #Rca1, Rpfc = (Munge.spiking.torate(@subset(spikes,:area .== ar), beh)
    #                for ar in ("CA1","PFC"))

    # Basic params
    # ----------------
    global embedding_overall
    using Munge.manifold
    load_manis(Main; feature_engineer, filt, distance, tag="$(animal)$(day).$(N)seg")
    embedding_overall = merge(embedding_overall, embedding)
end

embedding = embedding_overall


trans = :Matrix
transform(em) = if trans == :DimArray
    (em=em';DimArray(em, (Dim{:time}(beh.time[1:size(em,1)]), Dim{:comp}(1:size(em,2)))))
elseif trans == :Dataset
    (em=em';Dataset(eachcol(em)...))
else
    em'
end


# Which core would you like to work on?
#min_dist, n_neighbors, metric, dim = [0.3], [5,150], [:CityBlock, :Euclidean], 3
min_dist, n_neighbors, metric, dim = [0.3], [150], [:CityBlock], 3
K = filter(k->k.min_dist ∈ min_dist && k.n_neighbors ∈ n_neighbors && k.metric ∈ metric && k.dim == dim, keys(embedding))
embedding = Dict(k=>transform(embedding[k]) for k in K)

#embedding = 
