#!/bin/sh
#=
export LD_LIBRARY_PATH=/home/ryoung/miniconda3/envs/conda_jl/lib/
exec julia -J "/home/ryoung/Code/projects/goal-code/GFA-dependencies-sysimage.so" --project="/home/ryoung/Projects/goal-code/" "$0" -- $@
=#
using GoalFetchAnalysis, 
      GoalFetchAnalysis.Munge.manifold, 
      DrWatson, Revise
using GoalFetchAnalysis.Munge.manifold
using DataFramesMeta
using PyCall
@pyimport joblib as jl
@pyimport numpy as np

opt = Munge.manifold.parse(Main)
filt             = opt["filt"]
feature_engineer = opt["feature_engineer"]
distance         = opt["distance"]
sps              = opt["sps"]
splits           = opt["splits"]

animal, day = "RY22", 21
opt["animal"], opt["day"] = animal, day

for (animal, day) in DI.animal_set
    opt["animal"], opt["day"] = animal, day
    println("Running $animal $day")
    # FINISH loading libraries
    begin
        using Plots: StatsBase
        using Infiltrator, Serialization, Plots, ProgressMeter, PyCall, Distributed,
              ProgressMeter, ThreadSafeDicts, DataFramesMeta, Distances, StatsBase,
              SoftGlobalScope, Infiltrator, DimensionalData,  DataFramesMeta
        using DataStructures: OrderedDict
        import DIutils.namedtup: ntopt_string
        using PyCall
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
    sort!(spikes, :unit)
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
        behTrain = sort(behTrain, :time)
        spikesTrain = sort(spikesTrain, :time)
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
                 Munge.spiking.tocount(@subset(spikes,:area .== ar), beh)
                        for ar in ("CA1","PFC")
        )
        R = merge(R, 
                  Dict(
                      Symbol(lowercase(ar) * "_" * String(ct)) => 
                      (sub=@subset(spikes, :area .== ar, :celltype .== ct);
                       if !isempty(sub)
                           Munge.spiking.tocount(sub, beh)
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

    # Create a communication subspace channel
    ca1 = permutedims(Matrix(R[:ca1])[:,:,DIutils.na], [1,2,3])
    pfc = permutedims(Matrix(R[:pfc])[:,:,DIutils.na], [1,3,2])
    ca1pfc = ca1 .* pfc
    R = Dict{Symbol,DimArray}(R)
    R[:ca1pfc ] = DimArray(
        reshape(ca1pfc, size(ca1pfc,1),:), (:time, :ca1pfc))
    # Behavior variables to save
    props = [:time, :x, :y, :cuemem, :ha, :startWell, :stopWell, :traj, :epoch]
    savepath = datadir("exp_raw", "cebra")
    mkpath(savepath)
    (k,v) = first(R)

    for (k,v) in R
        k = string(k)
        println("Saving $k")
        B = Float64.(replace(Matrix(beh[!,props]), missing=>NaN))
        v = Matrix(v)
        sp = joinpath(savepath, "$animal-$day" * string(k)) * ".jl"
        py"""
        import joblib as jl
        import numpy as np
        jl.dump(dict(spikes=np.array($v, dtype=float), 
                position=np.array($B, dtype=float),
                animal=$animal, day=$day,dataset=$k),
                $sp,
                protocol=3)
        """
    end

end
