#!/bin/sh
#=
export LD_LIBRARY_PATH=/home/ryoung/miniconda3/envs/conda_jl/lib/
exec julia -J "/home/ryoung/Code/projects/goal-code/GFA-dependencies-sysimage.so" --project="/home/ryoung/Projects/goal-code/" "$0" -- $@
=#
using GoalFetchAnalysis, Munge.manifold, DrWatson, Revise
opt = Munge.manifold.parse()

# FINISH loading libraries
begin
    using Plots: StatsBase
    using Infiltrator, Serialization, Plots, ProgressMeter, PyCall, Distributed,
          ProgressMeter, ThreadSafeDicts, DataFramesMeta, Distances, StatsBase,
          SoftGlobalScope, Infiltrator, DimensionalData,  DataFramesMeta
    using DataStructures: OrderedDict
    import Utils.namedtup: ntopt_string
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
global areas            = (:ca1,:pfc)

# Load data
# ----------------


println("Loading")
@time global spikes, beh, ripples, cells = Load.load(animal, day)

println("Firing rate matrices")
R = Dict(Symbol(lowercase(ar))=>Munge.spiking.torate(@subset(spikes,:area .== ar), beh)
                for ar in ("CA1","PFC"))
zscoredimarray(x) = DimArray(hcat(zscore.(eachcol(x))...), x.dims)
R = merge(R,Dict(Symbol("z"*string(k))=>zscoredimarray(v) for (k,v) in R))


# Basic params
# ----------------

# Filter
println("Filtration?")
if filt !== nothing
    global filtstr = "filt=$filt"
    filters = Filt.get_filters()[filt]
    global beh,spikes = Utils.filtreg.filterAndRegister(beh, spikes; filters, 
    filter_skipmissingcols=true)
else
    global filtstr = "filt=nothing"
end
global festr   = feature_engineer === nothing ? "feature=nothing" : "feature=$feature_engineer"
global diststr = distance === nothing ? "distance=euclidean" : lowercase("distance=$distance")
@info "run info" filtstr festr diststr 

# Get sample runs
# ----------------
println("Generate partitions")
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
global N = splits * sps

global tag 
tag = "$(animal)$(day).$(N)seg"
println(tag)

# Get embeddings
# ----------------

#metrics      = unique((:CityBlock, :Euclidean,:Correlation,:Cosine))
#dimset       = (2,   3)
#min_dists    = (0.05,0.15,0.3)
#n_neighborss = (5,50,150,400)
min_dists, n_neighborss, metrics, dimset, features = [0.3], [5,150], [:CityBlock], 
                                                     [2,3], [:zscore]
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

trained_umap = em = sc = nothing
global steps, total = 0, (length(params) * length(datasets))

try
    for (metric,min_dist, n_neighbors, feature) in params
        for (area,dim,s) in datasets

            gc.collect()
            GC.gc(false)
            gc.collect()
            GC.gc(false)

            global steps += 1

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
            #pydecref(em)
            #pydecref(sc)
            pydecref(trained_umap)
            pydecref(fitter)
            fitter = trained_umap = em = sc = nothing
            next!(prog)
        end
    end

catch exception
    print(exception)
finally

    println("Finally section")
    println("Finally section")
    @info "Quit, probably GPU issue" animal day steps total steps/total
    # Store them for later
    using Munge.manifold
    savefile = path_manis(;filt,feature_engineer,tag)
    @info "save info" filt festr diststr savefile
    save_manis(;embedding, scores, inds_of_t, filt, feature_engineer, use_cuda, tag, splits, sps, N)
    #exit()
#end
end
