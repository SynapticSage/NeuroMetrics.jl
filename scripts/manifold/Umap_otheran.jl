quickactivate(expanduser("~/Projects/goal-code"))
using Infiltrator
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
using SoftGlobalScope, Infiltrator
using DimensionalData
use_cuda = true
if use_cuda
    using PyCall
    cuml = pyimport("cuml")
    @pyimport gc
    cuUMAP = cuml.manifold.umap.UMAP
end
using Munge.manifold

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

# DI data
# ----------------
try
animals = (("RY22", 21), ("RY16", 36))
#for (animal, day) in datasets 
(animal,day) = animals[1]

    @time global spikes, beh, ripples, cells = DI.load(animal, day)
    R = Dict(Symbol(lowercase(ar))=>Munge.spiking.torate(@subset(spikes,:area .== ar), beh)
                    for ar in ("CA1","PFC"))
    zscoredimarray(x) = DimArray(hcat(zscore.(eachcol(x))...), x.dims)
    R = merge(R,Dict(Symbol("z"*string(k))=>zscoredimarray(v) for (k,v) in R))


    # Basic params
    # ----------------
    global filt             = :all
    global areas            = (:ca1,:pfc)
    #distance        = :Mahalanobis
    global distance         = :many
    global feature_engineer = :many

    # Filter
    if filt !== nothing
        filtstr = "filt=$filt"
        filters = Filt.get_filters()[filt]
        global beh,spikes = Utils.filtreg.filterAndRegister(beh, spikes; filters, 
        filter_skipmissingcols=true)
    else
        filtstr = "filt=nothing"
    end
    global festr   = feature_engineer === nothing ? "feature=nothing" : "feature=$feature_engineer"
    global diststr = distance === nothing ? "distance=euclidean" : lowercase("distance=$distance")
    @info "run info" filtstr festr diststr 

    # Get sample runs
    # ----------------
    global splits, sampspersplit = 10, 10
    nsamp = Int(round(size(beh,1)/splits));
    δi    = Int(round(nsamp/sampspersplit))
    global inds_of_t = []
    for (split, samp) in Iterators.product(1:splits, 1:sampspersplit)
        start = (split-1) * nsamp + (samp-1) * δi + 1
        stop  = (split)   * nsamp + (samp-1) * δi + 1
        stop  = min(stop, size(beh,1))
        push!(inds_of_t, start:stop)
    end
    @info "coverage" nsamp/size(beh,1)*100
    global N = splits * sampspersplit
    global tag = "$(animal)$(day).$(N)seg"

    # Get embeddings
    # ----------------

    #metrics      = unique((:CityBlock, :Euclidean,:Correlation,:Cosine))
    #dimset       = (2,   3)
    #min_dists    = (0.05,0.15,0.3)
    #n_neighborss = (5,50,150,400)
    min_dists, n_neighborss, metrics, dimset, features = [0.3], [5,150], [:CityBlock], 
                                                         [2,3], [:raw,:zscore]
    #embedding,scores = Dict(), Dict()
    global embedding, scores = if isfile(path_manis(;filt,feature_engineer,tag))
        @info "loading prev data"
        data=load_manis(Main;filt,feature_engineer,tag);
        embedding, scores = data.embedding, data.scores
    else
        embedding, scores = Dict(), Dict()
    end
    @info "length of dicts" length(embedding) length(scores)
    using Infiltrator

    params   = collect(Iterators.product(metrics,min_dists, n_neighborss,features))
    datasets = collect(Iterators.product(areas, dimset, 1:length(inds_of_t)))
    prog = Progress(prod(length.((params,datasets))); desc="creating embeddings")

    #(min_dist, n_neighbors) = first(d

    trained_umap = em = sc = nothing
    for (metric,min_dist, n_neighbors, feature) in params
        for (area,dim,s) in datasets

            # Pre - process : Make key, do we process?
            key = (;area,dim,s,min_dist,n_neighbors,metric,feature)
            if key ∈ keys(embedding) && !(embedding[key] isa Future)
                @info "skipping" key
                next!(prog)
                continue
            else
                @info key
            end
            area = feature == :zscore ? Symbol("z"*string(area)) : area

            # Pre - process : Do we obtain a distance function?
            @debug "Getting distance metric"
            metric_str = if distance == :Mahalanobis
                @info "transforming Mahalanobis"
                Q = R[area]' * R[area]
                dist_func = getproperty(Distances, distance)(Matrix(Q[area]))
            elseif !use_cuda
                @info "transforming $distance"
                dist_func = getproperty(Distances, distance)
            else
                lowercase(string(metric))
            end


            input = Matrix(R[area]'[:,inds_of_t[s]])

            @debug "Processing"
            @time em, sc = if use_cuda # 1000x faster
                @debug "fitter"
                fitter=cuUMAP(n_neighbors=n_neighbors, min_dist=min_dist, 
                              n_components=3, metric=metric_str, 
                              target_metric="euclidean")
                @debug "fit"
                trained_umap = fitter.fit(input')
                @debug "transform"
                em = trained_umap.transform(input')
                @debug "trust"
                sc = cuml.metrics.trustworthiness(input', em, 
                                                  n_neighbors=n_neighbors)
                #@infiltrate
                em, sc
            else # julia is suprisingly slow here
                em = umap(input, dim; min_dist, n_neighbors, metric_str)
                skmani = pyimport("sklearn.manifold")
                @time sc = skmani.trustworthiness(input', em, 
                                               n_neighbors=n_neighbors, 
                                               metric=metric_str)
                em, sc
            end

            embedding[key] = em'
            scores[key] = sc
            trained_umap = em = sc = nothing
            gc.collect()
            next!(prog)
        end
    end

finally
    # Store them for later
    using Munge.manifold
    savefile = path_manis(;filt,feature_engineer,tag)
    @info "save info" filt festr diststr savefile
    save_manis(;embedding, scores, inds_of_t, filt, feature_engineer, use_cuda, tag, splits, sampspersplit, N)

    #exit()
#end
end

#data=load_manis(Main;filt,feature_engineer,tag);


#PL=@df sc scatter(:n_neighbors, :min_dist, :value;
#                  label="",xlabel="neighbors",ylabel="mindist",xscale=:log10)
#PL_nn=@df sc scatter(:n_neighbors, :value;
#                  label="",xlabel="neighbors",ylabel="mindist",xscale=:log10)
#PL_md=@df sc scatter(:min_dist, :value;
#                  label="",xlabel="neighbors",ylabel="mindist",xscale=:log10)
