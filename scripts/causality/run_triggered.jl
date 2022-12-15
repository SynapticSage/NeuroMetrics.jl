using DrWatson
using GoalFetchAnalysis
using Infiltrator
import Plot
using ThreadSafeDicts
include(scriptsdir("manifold","Umap_deserialize.jl"))

#  ================
# CAUSALITY AND MANFIOLD
#  ================

using Serialization, CausalityTools, Entropies
using Plots, DataFrames
using Statistics, NaNStatistics, HypothesisTests
using StatsPlots, JLD2
using Utils.namedtup: ntopt_string
using DataStructures: OrderedDict
using Munge.manifold, Munge.causal, Munge.triggering
using Utils.binning

include(scriptsdir("manifold","Umap_deserialize.jl"))

#  ================
# CAUSALITY AND MANFIOLD
#  ================


cor = Dict(0=>"correct", 1=>"error")
tsk = Dict(0=>"cue", 1=>"mem")
lab = OrderedDict([0,1]=>"CUE correct", [1,1]=>"MEM correct",
                  [0,0]=>"CUE error", [1,0]=>"MEM error")
conditionals = [:cuemem,:correct]

## ----------
## PARAMETERS
## ----------
animal, day = "RY16", 36

# Obtain params
# ---------------
@info "prop missinG"
esttype = :binned
est, params = get_est_preset(esttype)
params = (;params..., horizon=1:30, thread=false)
params = (;params..., binning=5, window=1.25)

# Obtain savefile
# ---------------
paramstr = Utils.namedtup.tostring(params)
N = 100
tagstr = if "tag" in propertynames(Main)
    "_$tag"
else
    "$animal$day.$(N)seg"
end


## ----------
## PARAMETERS
## ----------

"""
    obtain_triggered_causality(em, beh, props; grid_kws=nothing, 
        savefile,
        grd=(println("default binning compute");
             binning.get_grid(beh, props; grid_kws...)),
        checkpoint=ThreadSafeDict{NamedTuple, Bool}())

Gets the triggered causality, triggered by samples entering a bin. 
"""
function obtain_triggered_causality(em, beh, props; grid_kws=nothing, 
        savefile,
        grd=(println("default binning compute");
             binning.get_grid(beh, props; grid_kws...)),
        checkpoint=ThreadSafeDict{NamedTuple, Bool}())

    @info "obtain_local_binned_measure()" props savefile params params.window
    ## ------------------------------
    ## Setup frames and zone triggers
    ## ------------------------------
    #propstime = ["time",props...]
    EFF = EmbeddingFrameFetch(em, :area, beh, props; 
                              ordering=Dict(:area=>[:ca1,:pfc]))

    ## ------------------------
    ## Setup stores for results
    ## ------------------------
    if !isfile(savefile)
        ca1pfc = fill(NaN32,size(grd)..., length(EFF), params[:horizon].stop)
        pfcca1 = fill(NaN32,size(grd)..., length(EFF), params[:horizon].stop)
        ca1pfcσ² = fill(NaN32,size(grd)..., length(EFF), params[:horizon].stop)
        pfcca1σ² = fill(NaN32,size(grd)..., length(EFF), params[:horizon].stop)
        counts   = fill(Int16,size(grd)..., length(EFF), params[:horizon].stop)
    else
        @info "obtain_local_binned_measure, found save file, loading" savefile
        storage    = JLD2.jldopen(savefile, "r")
        checkpoint = merge(storage["checkpoint"], checkpoint)
        ca1pfc = storage["ca1pfc"] 
        pfcca1 = storage["pfcca1"] 
        ca1pfcσ² = storage["ca1pfcσ²"]
        pfcca1σ² = storage["pfcca1σ²"]
        counts = storage["counts"]
        close(storage)
    end
    

    ## ------------------------
    ## Compute
    ## ------------------------
    #(e,eff) = first(enumerate(EFF))

    @info "How much ground to cover" length(EFF) length(first(EFF))

    @showprogress "frames" for (e,eff) in enumerate(EFF[1:length(EFF)])
        @info e 
        zones = triggering.get_triggergen(params.window, grd, eff...)
        #(idx,data) = zones[50]
        Threads.@threads for (idx,data) in zones
            if isempty(data); continue; end
            checkpointkey = (;idx, e)
            if checkpointkey ∈ keys(checkpoint) && checkpoint[checkpointkey]
                continue
            else
                checkpoint[checkpointkey] = false
            end
            count = 0
            cpv = view(ca1pfc, idx..., e, :)
            pcv = view(pfcca1, idx..., e, :)
            cpvσ² = view(ca1pfcσ², idx..., e, :)
            pcvσ² = view(pfcca1σ², idx..., e, :)
            for D in data
                if isempty(D); continue; end
                ca1, pfc = map(d->transpose(hcat(d.data...)), D)
                try
                    cpt = causal.predictiveasymmetry(ca1, pfc; params...)
                    pct = causal.predictiveasymmetry(pfc, ca1; params...)
                    cpt,pct = convert(Vector{Float32}, cpt),
                              convert(Vector{Float32}, pct)
                    if all(isnan.(cpv)) && all(isnan.(pcv))
                        cpv .= 0
                        pcv .= 0
                    end
                    cpv .+= cpt
                    pcv .+= pct
                    cpvσ² .+= cpt.^2
                    pcvσ² .+= pct.^2
                catch err
                    showerror(stdout, err)
                end
                count += 1
            end
            ca1pfc[idx..., e, :] ./= count
            pfcca1[idx..., e, :] ./= count
            pfcca1σ²[idx..., e, :] ./= count
            ca1pfcσ²[idx..., e, :] ./= count
            counts[idx..., e, :] = count
            checkpoint[checkpointkey] = true
        end
        storage    = jldopen(savefile, "a")
        for field in ["ca1pfc","pfcca1","checkpoint","counts",
                      "ca1pfcσ²","pfcca1σ²"]
            delete!(storage,field)
            storage[field] = eval(Symbol(field))
        end
        close(storage)
    end


    # Better weighting
    weighting_trig = zeros(size(grd)..., length(EFF))
    times = []
    for (e,eff) in enumerate(EFF)
        zones = triggering.get_triggergen(1, grd, eff...)
        for (idx,data) in zones
            weighting_trig[idx..., e] += length(data)
            [append!(times, d[1].time) for d in data]
            #@infiltrate !isempty(data)
        end
        # @infiltrate all(weighting_trig[:,:,e] .== 0)
    end
    storage    = jldopen(savefile, "a")
    delete!(storage,"weighting_trig")
    storage["weighting_trig"] = weighting_trig
    (;ca1pfc, pfcca1, weighting_trig, checkpoint)
end
get_savefile(props;params=params) = 
    datadir("manifold","causal",
            "local_grid_cause_props=$(join(props,","))_$(Utils.namedtup.tostring(pop!(params,:thread)))_$tagstr.jld2")

# =====
# X - Y
# =====
grid_kws=(;widths=[5f0,5f0])
props = ["x","y"]
savefile = get_savefile(props)
ca1pfc, pfcca1, weighting_trig, checkpoint = 
        obtain_triggered_causality(em, beh, props; grid_kws, savefile)

JLD2.jldsave(savefile; params, est, grd, grid_kws, weighting_trig,
             props, ca1pfc, pfcca1, checkpoint)

# ==========================
# CUEMEM - CORRECT - HATRAJ
# ==========================
props = ["cuemem", "correct", "hatrajnum"]
grid_kws=(;widths=[1f0,1f0,1f0], 
           radiusinc=[0f0,0f0,0f0], 
           maxrad=[0.5f0,0.5f0, 0.5f0],
           radiidefault=[0.4f0,0.4f0,0.4f0],
           steplimit=1,
          )
grd = Utils.binning.get_grid(beh, props; grid_kws...)
savefile = get_savefile(props)

ca1pfc, pfcca1, weighting_trig, checkpoint = 
        obtain_triggered_causality(em, beh, props; grd, savefile)

hatrajcodex = Dict(r.hatrajnum=>r.hatraj for 
     r in eachrow(unique(beh[!,[:hatraj, :hatrajnum]])))
JLD2.jldsave(savefile; params, est, grd, grid_kws, weighting_trig,
             props, ca1pfc, pfcca1, checkpoint)

# ================================================
# CUEMEM - CORRECT - HATRAJ - startWell - stopWell
# ================================================
props = ["cuemem", "correct", "hatrajnum","startWell","stopWell"]
grid_kws=(;
           widths       = [1f0,1f0,1f0,1f0,1f0],
           radiusinc    = [0f0,0f0,0f0,0f0,0f0],
           maxrad       = [0.5f0,0.5f0, 0.5f0, 0.5f0, 0.5f0],
           radiidefault = [0.4f0,0.4f0,0.4f0,0.4f0,0.4f0],
           steplimit    = 1,
          )
grd = Utils.binning.get_grid(beh, props; grid_kws...)
savefile = get_savefile(props)

ca1pfc, pfcca1, weighting_trig, checkpoint = 
        obtain_triggered_causality(em, beh, props; grd, savefile)

hatrajcodex = Dict(r.hatrajnum=>r.hatraj for 
     r in eachrow(unique(beh[!,[:hatraj, :hatrajnum]])))
JLD2.jldsave(savefile; params, est, grd, grid_kws, weighting_trig,
             props, ca1pfc, pfcca1, checkpoint, hatrajcodex)

# =============================================================
# X - Y - cueme -  correct - hatrajnum - origin - destination 
# =============================================================
props = ["x", "y", "cuemem", "correct", "hatrajnum", "startWell", "stopWell"]
grid_kws =
        (;widths    = [4f0,4f0,1f0,1f0,1f0,1f0,1f0],
          radiusinc = [0.2f0,0.2f0,0f0,0f0,0f0,0f0,0f0],
          maxrad    = [6f0,6f0,0.4f0,0.4f0,0.4f0,0.4f0],
          radiidefault = [2f0,2f0,1f0,1f0,1f0,1f0,1f0]
         )
grd = Utils.binning.get_grid(beh, props; grid_kws...)
savefile = get_savefile(props)
ca1pfc, pfcca1, weighting_trig, checkpoint = 
        obtain_triggered_causality(em, beh, props; grd, savefile)

JLD2.jldsave(savefile; params, est, grd, grid_kws, weighting_trig,
             props, ca1pfc, pfcca1, checkpoint)

# ========================
# X - Y - START - STOP
# ========================
