using GoalFetchAnalysis, GoalFetchAnalysis.Munge.causal, 
GoalFetchAnalysis.Munge.manifold
using DIutils.namedtup: ntopt_string
import GoalFetchAnalysis.Plot
using DI.Labels
using DrWatson, Infiltrator, ProgressMeter, Serialization, CausalityTools,
Entropies, Plots, DataFrames, Statistics, NaNStatistics, HypothesisTests,
StatsPlots, JLD2, ThreadSafeDicts, SoftGlobalScope
using DataStructures: OrderedDict
function Base.fetch(X::AbstractDict)
    Dict(k=>(try;fetch(v);catch;missing;end) for (k,v) in X)
end

## ----------
## CONSTANTS
## ----------
corerr,tsk,cortsk = Labels.correct, Labels.cuemem, Labels.cortsk
## ----------
## PARAMS
## ----------
opt = Dict(
           :use_existing_predasym      => true, # use existing predasym[🔑]
           :skipexisting_predasym_keys => true, # skip existing predasym dict 
                #(if not, we can extend a key's results)
          )
PROPS = [[:cuemem, :correct],
         [:cuemem,:correct,:hatraj], 
         [:cuemem, :correct, :hatraj, :moving], 
         [:cuemem,:correct,:ha], 
         [:cuemem, :correct, :ha, :moving]]
datasets = (("RY22", 21, 6, :all), 
            ("RY16", 36, 6, :all))
#datasets = (("RY16", 36, 100, nothing),)
#datasets = (("RY22", 21, 100, nothing), )
@info "datasets" datasets
global predasym, savefile, loadfile = nothing, nothing, nothing
(animal, day, N, filt) = last(datasets)
(animal, day, N, filt) = first(datasets)


"""
    report_jobs_finished(C_pfcca1)
summarizes which threads are done, failed, or missing
# Arguments
- `C_pfcca1::Dict{Symbol,Dict{Symbol,Task}}`: a dictionary of dictionaries of
tasks
"""
function report_jobs_finished(C_pfcca1, 🔑, props, animal, day, N)
    itd(vv) = vv isa AbstractArray || istaskdone(vv)
    itf(vv) = vv isa AbstractArray || istaskfailed(vv)
    # Describe
    Idone=hcat([[!ismissing(vv) && itd(vv) for vv in values(v)]  
                for v in [V for V in values(C_pfcca1)]]...)
    Ifail=hcat([[!ismissing(vv) && itd(vv) for vv in values(v)] 
                for v in [V for V in values(C_pfcca1)]]...)
    Imissing=hcat([[!ismissing(vv) && itf(vv) for vv in values(v)] 
                   for v in [V for V in values(C_pfcca1)]]...)
    condition_data_plot = Plots.plot(Plots.heatmap(Idone, title="done"),
        heatmap(Ifail,title="fail"), heatmap(Imissing,title="missing"))
    setindex!.([info[🔑]], [Idone, Ifail, Imissing], ["done","fail","missing"])
    Plot.save("$animal.$day.$N.$props")
end


animal = day = N = filt = predasym = nothing
@showprogress "datasets" for (animal, day, N, filt) in datasets
    @info "dataset" animal day N filt
    predasym = nothing
    ## ----------
    ## PARAMETERS
    ## ----------
    areas = (:ca1,:pfc)
    distance = :many
    feature_engineer = :many # many | nothing
    esttype = :binned
    est, params = get_est_preset(esttype, horizon=1:60, thread=true, binning=7,
        window=1.25)
    # Loads
    manifold.load_manis_workspace(Main, animal, day; filt, 
                                  areas, distance, feature_engineer, 
                                  N)
    spikes, beh, ripples, cells  = DI.load(animal, day)
    # Save and load files
    savefile = get_alltimes_savefile(animal, day, N; params)
    loadfile = get_alltimes_savefile(animal, day, N; params,
        allowlowerhorizon=true)

    ## -----------------------------
    ## COMPUTE
    ## -----------------------------
    # GET PAIRED EMBEDDINGS
    # (each set of embeddings takes about 15 seconds)
    if isfile(loadfile)
        storage = JLD2.jldopen(loadfile, "r")
        CA1PFC, PFCCA1 = storage["CA1PFC"], storage["PFCCA1"]
        close(storage)
    else
        @assert !isempty(embedding)
        @time CA1PFC=Munge.causal.get_paired_embeddings(embedding, :zca1, :zpfc)
        @time PFCCA1=Munge.causal.get_paired_embeddings(embedding, :zpfc, :zca1)
        @time CA1CA1=Munge.causal.get_paired_embeddings(embedding, :zca1, :zca1)
        @time PFCPFC=Munge.causal.get_paired_embeddings(embedding, :zpfc, :zpfc)
        JLD2.jldsave(savefile; CA1PFC, PFCCA1, CA1CA1, PFCPFC, animal, day, N,
            params)
    end
    if CA1PFC == PFCCA1
        while (true); sleep(1); end
    end
    @assert first(values(CA1PFC)) != first(values(PFCCA1))
    # SETUP
    begin
        # Threading
        if params[:thread]
            Threads.nthreads() = 16
            print("Setting threads => ", Threads.nthreads())
            Dtype = ThreadSafeDict
        else
            Dtype = Dict
        end
        if opt[:use_existing_predasym]
            storagesave = isfile(savefile) ? jldopen(savefile, "r") : nothing
            storageload = isfile(loadfile) ? jldopen(loadfile, "r") : nothing
            if storagesave !== nothing && "predasym" in keys(storagesave)
                @info "loading prev predasym from storageSAVE"
                predasym = deepcopy(storagesave["predasym"])
            elseif storageload !== nothing && "predasym" in keys(storageload)
                @info "loading prev predasym from storageLOAD"
                predasym = deepcopy(storageload["predasym"])
            else
                predasym = nothing
            end
            storagesave === nothing ? nothing : close(storagesave)
            storageload === nothing ? nothing : close(storageload)
        end
        if predasym === nothing
            predasym = Dtype()
            predasym["alltimes"] = Dtype()
            for key in ("ca1pfc","pfcca1","ca1ca1","pfcpfc")
                predasym["alltimes"][key] = Dtype()
            end
        end
        info = Dtype()
        GC.gc()
    end
    # Function that handles computing conditional causal
    # measures
    """
        runconditional(ca1pfc, pfcca1, conditionals)
    Run conditional causal measures for `conditionals` on `ca1pfc` and `pfcca1`.
    # Arguments
    - `ca1pfc`: CA1 -> PFC embedding
    - `pfcca1`: PFC -> CA1 embedding
    - `conditionals`: list of symbols to condition on
    # Returns
    - `predasym`: Dict of D
    """
    function runconditional(ca1pfc, pfcca1, conditionals;
         ca1ca1=nothing, pfcpfc=nothing)
        if conditionals == [:cuemem, :correct]
            groups=[[1,1],[1,0],[0,1],[0,0]]
        elseif conditionals == [:cuemem, :correct, :hatraj]
            groups=[[cuemem, correct, hatraj] 
                    for cuemem in 0:1, correct in 0:1, 
                    hatraj in skipmissing(unique(beh.hatraj))]        
        elseif conditionals == [:cuemem, :correct, :hatrajnum]
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
        if ca1ca1 !== nothing
            @info "runconditional" "ca1ca1" conditionals
            conditional_pred_asym(ca1ca1, CA1CA1, beh, conditionals; groups, 
                                  inds_of_t, params...)
        end
        if pfcpfc !== nothing
            @info "runconditional" "pfcpfc" conditionals
            conditional_pred_asym(pfcpfc, PFCPFC, beh, conditionals; groups, 
                                  inds_of_t, params...)
        end
    end

    props, i = first(PROPS), 1
    props, i = last(PROPS),  2
    using SoftGlobalScope
    prog = Progress(length(PROPS), desc="Props")

    # Obtain conditional runs
    C_ca1pfc = C_pfcca1 = C_ca1ca1 = C_pfcpfc = nothing
    @softscope for (i,props) in enumerate(PROPS)
        🔑 = "props=" * replace(string(props),":"=>""," "=>"")
        @info "props" 🔑
        if 🔑 ∉ keys(predasym)
            predasym[🔑], info[🔑] = Dict{String,Any}(), Dict{String,Any}()
        elseif opt[:skipexisting_predasym_keys]
            print("🔑 in predasym, skipping...")
            continue
        end
        # Setup conditional storage
        if "ca1pfc" ∈ keys(predasym[🔑])
            C_ca1pfc = predasym[🔑]["ca1pfc"] 
            C_pfcca1 = predasym[🔑]["pfcca1"]
            C_ca1ca1 = predasym[🔑]["ca1ca1"]
            C_pfcpfc = predasym[🔑]["pfcpfc"]
        else
            C_ca1pfc = predasym[🔑]["ca1pfc"] = Dtype() 
            C_pfcca1 = predasym[🔑]["pfcca1"] = Dtype()
            C_ca1ca1 = predasym[🔑]["ca1ca1"] = Dtype()
            C_pfcpfc = predasym[🔑]["pfcpfc"] = Dtype()
        end
        # Run conditional
        runconditional(C_ca1pfc, C_pfcca1, props; 
                       ca1ca1=C_ca1ca1, pfcpfc=C_pfcpfc)
        run(`pushover-cli finished $animal iteration $i, key count = $(length(keys(C_ca1pfc)))`)
        # Fetch the results and store in master dict
        C_ca1pfc = fetch(C_ca1pfc)
        C_pfcca1 = fetch(C_pfcca1)
        C_ca1ca1 = fetch(C_ca1ca1)
        C_pfcpfc = fetch(C_pfcpfc)
        predasym[🔑]["ca1pfc"] = C_ca1pfc
        predasym[🔑]["pfcca1"] = C_pfcca1
        predasym[🔑]["ca1ca1"] = C_ca1ca1
        predasym[🔑]["pfcpfc"] = C_pfcpfc
        predasym = fetch(predasym)
        report_jobs_finished(C_pfcca1, 🔑, props, animal, day, N)
        # ==========
        # Checkpoint
        # ==========
        using JLD2
        S=JLD2.jldsave(savefile; CA1PFC, PFCCA1, animal, day, N, params, est, 
                       predasym, info)
        run(`pushover-cli saved $animal iteration $i, key count = $(length(keys(C_ca1pfc)))`)
        next!(prog)
    end

    # COMPUTE ALL TIMES
    if isempty(predasym["alltimes"]["ca1pfc"]) || !opt[:skipexisting_predasym_keys]
        predictiveasymmetry!(predasym["alltimes"]["ca1pfc"], CA1PFC; params...)
        predictiveasymmetry!(predasym["alltimes"]["pfcca1"], PFCCA1; params...)
        predictiveasymmetry!(predasym["alltimes"]["ca1ca1"], CA1CA1; params...)
        predictiveasymmetry!(predasym["alltimes"]["pfcpfc"], PFCPFC; params...)
        predasym["alltimes"]["ca1pfc"] = fetch(predasym["alltimes"]["ca1pfc"])
        predasym["alltimes"]["pfcca1"] = fetch(predasym["alltimes"]["pfcca1"])
        predasym["alltimes"]["ca1ca1"] = fetch(predasym["alltimes"]["ca1ca1"])
        predasym["alltimes"]["pfcpfc"] = fetch(predasym["alltimes"]["pfcpfc"])
        run(`pushover-cli finished $animal alltimes`)
    end

    # Checkpoint
    using JLD2
    S=JLD2.jldsave(savefile; CA1PFC, PFCCA1, animal, day, N, params, est, 
                   predasym, info)

end # Animal dataset loop!
