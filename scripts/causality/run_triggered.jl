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

using Munge.causal
get_savefile(props;params=params) = 
    datadir("manifold","causal",
            "local_grid_cause_props=$(join(props,","))_$(Utils.namedtup.tostring(pop!(params,:thread)))_$tagstr.jld2")

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

Munge.causal.obtain_triggered_causality(em, beh, props, params; 
                                 grd, savefile, checkpointmod=5)
hatrajcodex = Dict(r.hatrajnum=>r.hatraj for 
     r in eachrow(unique(beh[!,[:hatraj, :hatrajnum]])))

storage = jldopen(savefile,"a");
storage["grd"]         = grd
storage["params"]      = params
storage["grid_kws"]    = grid_kws
storage["hatrajcodex"] = hatrajcodex
close(storage)

# =============================================================
# X - Y - cueme -  correct - hatrajnum - origin - destination 
# =============================================================
props = ["x", "y", "cuemem", "correct", "hatrajnum", "startWell", "stopWell"]
grid_kws =
        (;widths    = [4f0,4f0,1f0,1f0,1f0,1f0,1f0],
          radiusinc = [0.2f0,0.2f0,0f0,0f0,0f0,0f0,0f0],
          maxrad    = [6f0,6f0,0.4f0,0.4f0,0.4f0,0.4f0,0.4f0],
          radiidefault = [2f0,2f0,0.4f0,0.4f0,0.4f0,0.4f0,0.4f0],
          steplimit=3,
         )
savefile = get_savefile(props)
storage = jldopen(savefile,"a");
if "grd" ∈ keys(storage)
    grd = storage["grd"]
    checkpoint = storage["checkpoint"]
    @assert params == storage["params"]
else
    grd = Utils.binning.get_grid(beh, props; grid_kws...)
    storage["grd"] = grd
    storage["grid_kws"] = grid_kws
    storage["params"] = params
    storage["checkpoint"] = checkpoint = ThreadSafeDict{NamedTuple, Bool}()
end

close(storage)
Munge.causal.obtain_triggered_causality(em, beh, props, params; 
                                        grd, savefile, checkpointmod=5,
                                        floattype=Float16, inttype=Int32,
                                        nan=NaN16, checkpoint
                                       )

JLD2.jldsave(savefile; params, est, grd, grid_kws, weighting_trig,
             props, ca1pfc, pfcca1, checkpoint)


# ==========================
# CUEMEM - CORRECT - HATRAJ
# ==========================
props = ["cuemem", "correct"]
grid_kws=(;widths=[1f0,1f0], 
           radiusinc=[0f0,0f0], 
           maxrad=[0.5f0,0.5f0],
           radiidefault=[0.4f0,0.4f0],
           steplimit=1,
          )

grd = Utils.binning.get_grid(beh, props; grid_kws...)
savefile = get_savefile(props)
obtain_triggered_causality(em, beh, props; grd, savefile, checkpointmod=5)

hatrajcodex = Dict(r.hatrajnum=>r.hatraj for 
     r in eachrow(unique(beh[!,[:hatraj, :hatrajnum]])))
JLD2.jldsave(savefile; params, est, grd, grid_kws, weighting_trig,
             props, ca1pfc, pfcca1, checkpoint)
# ========================
# X - Y - START - STOP
# ========================
## =====
## X - Y
## =====
#grid_kws=(;widths=[5f0,5f0])
#props = ["x","y"]
#savefile = get_savefile(props)
#ca1pfc, pfcca1, weighting_trig, checkpoint = 
#        obtain_triggered_causality(em, beh, props; grid_kws, savefile)
#
#JLD2.jldsave(savefile; params, est, grd, grid_kws, weighting_trig,
#             props, ca1pfc, pfcca1, checkpoint)

