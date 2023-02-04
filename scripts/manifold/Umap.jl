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
beh.index = 1:size(beh,1)

# ----------------
# Basic params
# ----------------
global spikesTrain = spikesTest = spikes
global behTrain = behTest = beh

# Filter
println("Filtration?")
if filt !== nothing
    global filtstr = "filt=$filt"
    filters = Filt.get_filters()[filt]
    global behTrain, spikesTrain =
                        Utils.filtreg.filterAndRegister(copy(behTrain), copy(spikesTrain); filters,
                                                       filter_skipmissingcols=true)
else
    global filtstr = "filt=nothing"
end
global festr   = feature_engineer === nothing ? "feature=nothing" : "feature=$feature_engineer"
global diststr = distance === nothing ? "distance=euclidean" : lowercase("distance=$distance")
@info "run info" filtstr festr diststr 

# Firing Rates
# ------------
println("Firing rate matrices")
function get_R(beh, spikes)
    R = Dict(
             Symbol(lowercase(ar))=>Munge.spiking.torate(@subset(spikes,:area .== ar), beh)
                    for ar in ("CA1","PFC")
            )
    zscoredimarray(x) = DimArray(hcat(zscore.(eachcol(x))...), x.dims)
    merge(R,Dict(Symbol("z"*string(k))=>zscoredimarray(v) for (k,v) in R))
end
R      = get_R(beh, spikes)
# Rtrain = get_R(behTrain, spikesTrain)
Rtrain_inds = Dict(k=>Utils.searchsortednearest.([beh.time], behTrain.time) for (k,_) in R)
Rtrain = Dict(
              k=>v[inds,:] for ((k,v),(k,inds)) in zip(R,Rtrain_inds)
             )

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
datasets = collect(Iterators.product(areas, dimset, 1:length(inds_of_t)))
prog = Progress(prod(length.((params,datasets))); desc="creating embeddings")
global steps, total = 0, (length(params) * length(datasets))

trained_umap = em = sc = nothing
exception_triggered = false

try

    for (metric,min_dist, n_neighbors, feature) in params
        for (rate_data,dim,s) in datasets

            gc.collect()
            GC.gc(false)
            gc.collect()
            GC.gc(false)

            global steps += 1

            # Pre - process : Make key, do we process?
            key = (;rate_data,dim,s,min_dist,n_neighbors,metric,feature)
            if key ∈ keys(embedding) && !(embedding[key] isa Future)
                @info "skipping" key
                next!(prog)
                continue
            else
                @info key
            end
            
            area = feature == :zscore ? Symbol("z"*string(rate_data)) : rate_data;

            # Pre - process : Do we obtain a distance function?
            @debug "Getting distance metric"
            metric_str = if distance == :Mahalanobis
                @info "transforming Mahalanobis"
                Q = Rtrain[area]' * Rtrain[area];
                dist_func = getproperty(Distances, distance)(Matrix(Q[area]));
            elseif !use_cuda
                @info "transforming $distance"
                dist_func = getproperty(Distances, distance)
            else
                lowercase(string(metric))
            end

            train = Matrix(Rtrain[area]'); # train on everything with filters
            #train = Matrix(Rtrain[area]'[:, # train on just this set of indices (dynamic manifold)
                                         # Rtrain_matching[area, inds_of_t[s]]
                                        # ]); 
            input = Matrix(R[area]'[:, # test on this set minus filters
                                    inds_of_t[s]
                                   ]);
            # input = Matrix(R[area]'); # test on everything minus filters

            @debug "Processing"
            @time em, sc = if use_cuda # 1000x faster
                @debug "fitter"

                dim = 3
                fitter=cuUMAP(n_neighbors=n_neighbors, min_dist=min_dist, 
                              n_components=dim, metric=metric_str, local_connectivity=dim,
                              target_metric="euclidean", n_epochs=1500);
                # local_connectivity: int (optional, default 1)
                #     The local connectivity required -- i.e. the number of nearest
                #     neighbors that should be assumed to be connected at a local level.
                #     The higher this value the more connected the manifold becomes
                #     locally. In practice this should be not more than the local intrinsic
                #     dimension of the manifold.
                # n_epochs larger = more accurate, 500 for small, 200 for large
                randsamp(x,n) = x[:,Random.randperm(size(x,2))[1:n]]
                @debug "fit"
                T = train
                trained_umap = fitter.fit(T');
                @debug "transform"
                # I = randsamp(input, 50_000)
                I = input
                em = trained_umap.transform(I');
                @debug "trust"
                sc = cuml.metrics.trustworthiness(I', em, n_neighbors=n_neighbors)
                # KEEP THIS HERE FOR INTRA-SCRIPT TESTING
                # k,w=0.5,20;
                # import Random
                # E = em[Random.randperm(size(em,1))[1:14_000],:];
                # begin
                #     s=scatter(eachcol(E)...; alpha=0.1 * k, markersize=2);
                #     plot!(eachcol(E)..., alpha=0.01 * k); ylims!(-w, w);
                #     xlims!(-w,w); zlims!(-w,w);
                #     s
                # end
                # t=@task Plot.stereoscopicgif( eachcol(E)... ;deltaangle=5, xlims=(-w,w), ylim=(-w,w), zlim=(-w,w), alphaasync =0.1*k)
                em, sc
            else # julia is suprisingly slow here
                em = umap(input, dim; min_dist, n_neighbors, metric_str)
                skmani = pyimport("sklearn.manifold")
                @time sc = skmani.trustworthiness(input', em, 
                                               n_neighbors=n_neighbors, 
                                               metric=metric_str)
                em, sc
            end;

            embedding[key] = em';
            scores[key] = sc;
            #pydecref(em)
            #pydecref(sc)
            pydecref(trained_umap);
            pydecref(fitter);
            fitter = trained_umap = em = sc = nothing;
            next!(prog)
        end
    end

catch exception
    exception_triggered = true
    print(exception)
finally

    println("Finally section")
    if exception_triggered
        @info "Quit, probably GPU issue" animal day steps total steps/total
    else
        @info "Finished! 😺" animal day steps total steps/total
    end
    # Store them for later
    using Munge.manifold
    savefile = path_manis(;filt,feature_engineer,tag)
    @info "save info" filt festr diststr savefile
    save_manis(;embedding, scores, inds_of_t, filt, feature_engineer, use_cuda, tag, splits, sps, N)
    #exit()
#end
end
