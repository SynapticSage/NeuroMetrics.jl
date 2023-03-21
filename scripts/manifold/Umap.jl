#!/bin/sh
#=
export LD_LIBRARY_PATH=/home/ryoung/miniconda3/envs/conda_jl/lib/
exec julia -J "/home/ryoung/Code/projects/goal-code/GFA-dependencies-sysimage.so" --project="/home/ryoung/Projects/goal-code/" "$0" -- $@
=#
using GoalFetchAnalysis, 
      GoalFetchAnalysis.Munge.manifold, 
      DrWatson, Revise
using GoalFetchAnalysis.Munge.manifold
opt = Munge.manifold.parse()
filt             = opt["filt"]
feature_engineer = opt["feature_engineer"]
distance         = opt["distance"]
sps              = opt["sps"]
splits           = opt["splits"]

# Sets to explore
# ---------------
#metrics      = unique((:CityBlock, :Euclidean,:Correlation,:Cosine))
#dimset       = (2,   3)
#min_dists    = (0.05,0.15,0.3)
# #n_neighborss = (5,50,150,400)
# min_dists, n_neighborss, metrics, dimset, features = [0.5], [400],
# [:Euclidean], [2,3], [:zscore]
# Input Wenbo's settings
min_dists, n_neighborss, metrics, dimset, features = [0.3], [100], 
                                        [:Cosine], [6], [:zscore]
negative_sample_rate = 100

function keyfunc(;dataset,dim,s,min_dist,n_neighbors,metric,feature)
    (;dataset,dim,s,min_dist,n_neighbors,metric,feature)
end

# FINISH loading libraries
begin
    using Plots: StatsBase
    using Infiltrator, Serialization, Plots, ProgressMeter, PyCall, Distributed,
          ProgressMeter, ThreadSafeDicts, DataFramesMeta, Distances, StatsBase,
          SoftGlobalScope, Infiltrator, DimensionalData,  DataFramesMeta
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

# ----------------
# DI data
# ----------------

println("Loading")
@time global spikes, beh, ripples, cells = DI.load(opt["animal"], opt["day"])
cells, spikes = DIutils.filtreg.register(cells, spikes; on="unit",
    transfer=["celltype"])
beh.index = 1:size(beh,1)

# ----------------
# Basic params
# ----------------
global spikesTrain = spikesTest = spikes
global behTrain = behTest = beh

# Filter
println("Filtration?")
if opt["filt"] !== nothing
    global filtstr = "filt=$(opt["filt"])"
    filters = Filt.get_filters()[opt["filt"]]
    global behTrain, spikesTrain =
                        DIutils.filtreg.filterAndRegister(copy(behTrain),
            copy(spikesTrain); filters, filter_skipmissingcols=true)
else
    global filtstr = "filt=nothing"
end
global festr   = opt["feature_engineer"] === nothing ? "feature=nothing" :
    "feature=$(opt["feature_engineer"])"
global diststr = opt["distance"] === nothing ? "distance=euclidean" :
    lowercase("distance=$(opt["distance"])")
@info "run info" filtstr festr diststr 

# Firing Rates
# ------------
println("Firing rate matrices")

function get_R(beh, spikes; gaussian=0.25)
    R = Dict(
             Symbol(lowercase(ar)) =>
             Munge.spiking.torate(@subset(spikes,:area .== ar), beh; gaussian)
                    for ar in ("CA1","PFC")
    )
    R = merge(R, 
              Dict(
                  Symbol(lowercase(ar) * "_" * String(ct)) => 
                  (sub=@subset(spikes, :area .== ar, :celltype .== ct);
                   if !isempty(sub)
                       Munge.spiking.torate(sub, beh; gaussian=gaussian)
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
R      = get_R(beh, spikes)

# ----------------
#
# ----------------
# Rtrain = get_R(behTrain, spikesTrain)
Rtrain_inds = Dict(k=>DIutils.searchsortednearest.([beh.time], behTrain.time) for (k,_) in R)
Rtrain = Dict(
              k=>v[inds,:] for ((k,v),(k,inds)) in zip(R,Rtrain_inds)
             )

# ----------------
# Get sample runs
# ----------------
println("Generate partitions")
nsamp = Int(round(size(beh,1)/opt["splits"]));
δi    = Int(round(nsamp/opt["sps"]))
global inds_of_t = []
for (split, samp) in Iterators.product(1:opt["splits"], 1:opt["sps"])
    start = (split-1) * nsamp + (samp-1) * δi + 1
    stop  = (split)   * nsamp + (samp-1) * δi + 1
    stop  = min(stop, size(beh,1))
    push!(inds_of_t, start:stop)
end

println("Describ partitions")
@info "coverage" nsamp/size(beh,1)*100
global N = opt["splits"] * opt["sps"]

global tag = "$(opt["animal"])$(opt["day"]).$(N)seg"
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
fitter_params = nothing


function save_results()
    # Store them for later
    savefile = path_manis(;filt,feature_engineer,tag)
    @info "save info" filt festr diststr savefile
    try
        save_manis(;embedding, scores, inds_of_t, filt, feature_engineer, use_cuda,
                   tag, splits, sps, N, fitter_params)
        printstyled("FINISHED UMAP.JL"; blink=true, color=:green)
        cmd = 
        `pushover-cli "Finished UMAP for animal $animal day $day filter $filt"`
        run(cmd)
    catch
        printstyled("FAILED SAVE UMAP.JL"; blink=true, color=:light_red, reverse=true)
        cmd = 
        `pushover-cli "Unfinished UMAP for animal $animal day $day filter $filt"`
        run(cmd)
    finally
        keyset = [keyfunc(;metric,min_dist,n_neighbors, feature, dataset,dim,s)  for
         (metric,min_dist, n_neighbors, feature) in params, (dataset,dim,s) in datasets]
        failed_keys = [k for k in keyset if k ∉ keys(embedding)]
        succes_keys = [k for k in keyset if k ∈ keys(embedding)]
        @info "failed keys => $failed_keys"
    end
end

steps = 0;
try

    (metric,min_dist, n_neighbors, feature) = first(params)
    for (metric,min_dist, n_neighbors, feature) in params

        (dataset,dim,s) = first(datasets)
        for (dataset,dim,s) in datasets

            steps += 1

            gc.collect()
            GC.gc(false)
            gc.collect()
            GC.gc(false)

            global steps += 1
            key = keyfunc(;dataset,dim,s,min_dist,n_neighbors,metric,feature)
            # Pre - process : Make key, do we process?
            if key ∈ keys(embedding) && !(embedding[key] isa Future)
                @info "skipping" key
                next!(prog)
                continue
            else
                @info key
            end
            
            if feature == :zscore && !(startswith(string(dataset), "z"))
                continue
            end

            # Pre - process : Do we obtain a distance function?
            @debug "Getting distance metric"
            metric_str = if opt["distance"] == :Mahalanobis
                @info "transforming Mahalanobis"
                Q = Rtrain[dataset]' * Rtrain[dataset];
                dist_func = getproperty(Distances, opt["distance"])(Matrix(Q[dataset]));
            elseif !use_cuda
                @info "transforming $(opt["distance"])"
                dist_func = getproperty(Distances, opt["distance"])
            else
                lowercase(string(metric))
            end

            # train = Matrix(Rtrain[dataset]'); # train on everything with filters
            train = Matrix(Rtrain[dataset]'[:, # train on just this set of indices (dynamic manifold)
                                         Rtrain_matching[dataset, inds_of_t[s]]
                                        ]); 
            input = Matrix(R[dataset]'[:, # test on this set minus filters
                                    inds_of_t[s]
                                   ]);
            # input = Matrix(R[dataset]'); # test on everything minus filters

            @debug "Processing"
            @time em, sc = if use_cuda # 1000x faster
                @debug "fitter"

                T = train
                I = input
                # T = project_onto_behavior(T', beh, inds_of_t[s], beh_vars)'
                # I = project_onto_behavior(I, beh, inds_of_t[s], beh_vars)

                fitter=cuUMAP(n_neighbors=n_neighbors, min_dist=min_dist,
                n_components=dim, metric="cosine", local_connectivity=10,
                repulsion_strength=1, negative_sample_rate=negative_sample_rate,
                target_metric="euclidean", n_epochs=10_000);
                fitter_params = fitter.get_params()
                # local_connectivity: int (optional, default 1)
                # The local connectivity required -- i.e. the number of
                # nearest neighbors that should be assumed to be connected at a
                # local level. The higher this value the more connected the
                # manifold becomes locally. In practice this should be not more
                # than the local intrinsic dimension of the manifold.
                # n_epochs larger = more accurate, 500 for small, 200 for large
                @debug "fit"
                py"def fit_and_catch_errors(T):
                    try:
                        return $(fitter).fit(T)
                    except (ValueError, ArgumentError) as e:
                        print(e)
                        return None"
                fit_func = py"fit_and_catch_errors"
                # trained_umap = fitter.fit(T');
                trained_umap = fit_func(T')
                @debug "transform"
                # randsamp(x,n) = x[:,Random.randperm(size(x,2))[1:n]]
                # I = randsamp(input, 50_000)
                em = trained_umap.transform(I');
                @debug "trust"
                sc = cuml.metrics.trustworthiness(I', em, n_neighbors=n_neighbors)
                # t=@task Plot.stereoscopicgif( eachcol(E)... ;deltaangle=5, xlims=(-w,w), ylim=(-w,w), zlim=(-w,w), alphaasync =0.1*k)
                em, sc

            else # julia is suprisingly slow here
                import UMAP
                em = UMAP.umap(input, dim; min_dist, n_neighbors, metric_str)
                skmani = pyimport("sklearn.manifold")
                @time sc = skmani.trustworthiness(input', em, 
                                               n_neighbors=n_neighbors, 
                                               metric=metric_str)
                println("dataset $dataset, s $s, trustworthiness $sc")
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

            if steps % opt["save_frequency"] == 0
                save_results()
            end

        end
    end

catch exception
    exception_triggered = true
    print(exception)
    save_results()
finally
    println("Finally section")
    if exception_triggered
        @info "Quit, probably GPU issue" opt steps total steps/total
    else
        @info "Finished! 😺" opt steps total steps/total
    end
    save_results()
end

