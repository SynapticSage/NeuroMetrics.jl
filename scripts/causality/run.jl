using GoalFetchAnalysis
import Plot
using Utils.namedtup: ntopt_string
using Munge.causal, Munge.manifold

using DrWatson
using Infiltrator, ProgressMeter
using Serialization
using CausalityTools
using Entropies
using DataStructures: OrderedDict
using Plots
using DataFrames
using Statistics, NaNStatistics, HypothesisTests
using StatsPlots
using JLD2
using DiskBackedDicts, ThreadSafeDicts
using SoftGlobalScope
function Base.fetch(X::AbstractDict)
    Dict(k=>(try;fetch(v);catch;missing;end) for (k,v) in X)
end

## ----------
## CONSTANTS
## ----------
corerr,tsk,cortsk = Munge.behavior.cor, Munge.behavior.tsk, 
              Munge.behavior.cortsk


## ----------
## PARAMS
## ----------
PROPS = [[:cuemem, :correct],[:cuemem,:correct,:hatraj]]
datasets = (("RY22", 21, 100, nothing), ("RY16", 36, 100, nothing))
datasets = (("RY16", 36, 100, nothing),)
opts = Dict(
            :skipifproc => true
           )

global predasym = nothing
(animal, day, N, filt) = last(datasets)
(animal, day, N, filt) = first(datasets)


@showprogress "datasets" for (animal, day, N, filt) in datasets

    @info "dataset" animal day N filt

    ## ----------
    ## PARAMETERS
    ## ----------
    areas = (:ca1,:pfc)
    distance = :many
    feature_engineer = :many # many | nothing
    esttype = :binned
    est, params = get_est_preset(esttype)
    params = (;params..., horizon=1:30, thread=true)
    params = (;params..., binning=7, window=1.25)
    manifold.load_manis_workspace(Main, animal, day; filt, 
                                  areas, distance, feature_engineer, 
                                  N)
    spikes, beh, ripples, cells  = Load.load(animal, day)
    savefile = get_alltimes_savefile(animal, day, N; params)

    ## -----------------------------
    ## COMPUTE
    ## -----------------------------

    # GET PAIRED EMBEDDINGS
    # (each set of embeddings takes about 15 seconds)
    if isfile(savefile)
        storage = JLD2.jldopen(savefile, "r")
        CA1PFC, PFCCA1 = storage["CA1PFC"], storage["PFCCA1"]
        close(storage)
    else
        @assert !isempty(embedding)
        @time CA1PFC = Munge.causal.get_paired_embeddings(embedding, :ca1, :pfc)
        @time PFCCA1 = Munge.causal.get_paired_embeddings(embedding, :pfc, :ca1)
        JLD2.jldsave(savefile; CA1PFC, PFCCA1, animal, day, N, params)
    end

    # SETUP
    usecheckpoint = true
    begin
        # Threading
        if params[:thread]
            Threads.nthreads() = 16
            print("Threads => ", Threads.nthreads())
            Dtype = ThreadSafeDict
        else
            Dtype = Dict
        end
        storage = jldopen(savefile, "r")
        if usecheckpoint && "predasym" in keys(storage)
            predasym = storage["predasym"]
        else
            predasym = Dtype()
            predasym["alltimes"] = Dtype()
            for key in ("ca1pfc","pfcca1")
                predasym["alltimes"][key] = Dtype()
            end
        end
        close(storage)
        info = Dtype()
        GC.gc()
    end

    # Function that handles computing conditional causal
    # measures
    function runconditional(ca1pfc, pfcca1, conditionals)
        if conditionals == [:cuemem, :correct]
            groups=[[1,1],[1,0],[0,1],[0,0]]
        elseif conditionals == [:cuemem, :correct, :hatraj]
            groups=[[cuemem, correct, hatraj] 
                    for cuemem in 0:1, correct in 0:1, 
                    hatraj in skipmissing(unique(beh.hatraj))]        elseif conditionals == [:cuemem, :correct, :hatrajnum]
            groups=[[cuemem, correct, hatraj] 
                    for cuemem in 0:1, correct in 0:1, 
                    hatraj in skipmissing(unique(beh.hatrajnum))]
        else
            groups=nothing
        end
        @info "runconditional" "ca1pfc" conditionals
        conditional_pred_asym(ca1pfc, CA1PFC, beh, conditionals; groups, 
                              inds_of_t, params...)
        @info "runconditional" "pfcca1" conditionals
        conditional_pred_asym(pfcca1, PFCCA1, beh, conditionals; groups, 
                              inds_of_t, params...)
    end

    # Obtain conditional runs
    props, i = first(PROPS), 1
    props, i = last(PROPS), 2
    using SoftGlobalScope
    prog = Progress(length(PROPS), desc="Props")
    @softscope for (i,props) in enumerate(PROPS)

        🔑 = "props=" * replace(string(props),":"=>""," "=>"")
        @info "props" 🔑
        if 🔑 ∉ keys(predasym)
            predasym[🔑], info[🔑] = Dict{String,Any}(), Dict{String,Any}()
        else
            print("🔑 in predasym, skipping...")
            continue
        end
        @infiltrate

        C_ca1pfc = predasym[🔑]["ca1pfc"] = Dtype() 
        C_pfcca1 = predasym[🔑]["pfcca1"] = Dtype()

        runconditional(C_ca1pfc, C_pfcca1, props)
        run(`pushover-cli finished $animal iteration $i, key count = $(length(keys(C_ca1pfc)))`)

        # Describe
        Idone=hcat([[!ismissing(vv) && istaskdone(vv) for vv in values(v)]  
                    for v in [V for V in values(C_pfcca1)]]...)
        Ifail=hcat([[!ismissing(vv) && istaskfailed(vv) for vv in values(v)] 
                    for v in [V for V in values(C_pfcca1)]]...)
        Imissing=hcat([[!ismissing(vv) && istaskfailed(vv) for vv in values(v)] 
                       for v in [V for V in values(C_pfcca1)]]...)
        condition_data_plot = Plots.plot(Plots.heatmap(Idone, title="done"), 
                                   heatmap(Ifail,title="fail"), heatmap(Imissing,title="missing"))
        setindex!.([info[🔑]], [Idone, Ifail, Imissing], ["done","fail","missing"])
        Plot.save("$animal.$day.$N.$props")

        C_ca1pfc = fetch(C_ca1pfc)
        C_pfcca1 = fetch(C_pfcca1)
        predasym[🔑]["ca1pfc"] = C_ca1pfc
        predasym[🔑]["pfcca1"] = C_pfcca1
        predasym = fetch(predasym)

        # Checkpoint
        using JLD2
        S=JLD2.jldsave(savefile; CA1PFC, PFCCA1, animal, day, N, params, est, 
                       predasym, info)

        next!(prog)

    end

    # COMPUTE ALL TIMES
    begin
        predictiveasymmetry!(predasym["alltimes"]["ca1pfc"], CA1PFC; params...)
        predictiveasymmetry!(predasym["alltimes"]["pfcca1"], PFCCA1; params...)
        predasym["alltimes"]["ca1pfc"] = fetch(predasym["alltimes"]["ca1pfc"])
        predasym["alltimes"]["pfcca1"] = fetch(predasym["alltimes"]["pfcca1"])
    end

    # Checkpoint
    using JLD2
    S=JLD2.jldsave(savefile; CA1PFC, PFCCA1, animal, day, N, params, est, 
                   predasym, info)

end
