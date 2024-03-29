#  ==================================
# CAUSALITY AND MANFIOLD :: Triggered
#  ==================================
using DrWatson
using Infiltrator
using ThreadSafeDicts
using Serialization, CausalityTools, Entropies
using Plots, DataFrames
using Statistics, NaNStatistics, HypothesisTests
using StatsPlots, JLD2
using DataStructures: OrderedDict

using GoalFetchAnalysis
using Utils.namedtup: ntopt_string
import Plot
using Munge.manifold, Munge.causal, Munge.triggering
using Utils.binning
using Munge.causal



## ----------
## CONSTANTS
## ----------
cor,tsk,cortsk = Munge.behavior.cor, Munge.behavior.task, 
              Munge.behavior.cortsk

## ----------
## PARAMETERS
## ----------
animal, day, N, filt = "RY16", 36, 100, nothing
areas = (:ca1,:pfc)
distance = :many
feature_engineer = :many # many | nothing
esttype = :binned
est, params = get_est_preset(esttype)
params = (;params..., horizon=1:30, thread=false)
params = (;params..., binning=5, window=1.25)
manifold.load_manis_workspace(Main, animal, day; filt, 
                              areas, distance, feature_engineer, 
                              N)

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
savefile = get_trigger_savefile(props)

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

savefile = get_trigger_savefile(props)
storage = jldopen(savefile,"a");
if "grd" ∈ keys(storage)
    grd = storage["grd"]
    if "checkpoint" in keys(storage)
        checkpoint = storage["checkpoint"]
    else
        checkpoint = ThreadSafeDicts{NamedTuple, Bool}()
    end
    @assert params == storage["params"]
else
    grd = Utils.binning.get_grid(beh, props; grid_kws...)
    storage["grd"]      = grd
    storage["grid_kws"] = grid_kws
    storage["params"]   = params
    storage["checkpoint"] = checkpoint = ThreadSafeDict{NamedTuple, Bool}()
end

close(storage)
Munge.causal.obtain_triggered_causality_binavg(em, beh, props, params; 
                                        grd, savefile, checkpointmod=5,
                                        floattype=Float16, inttype=Int32,
                                        nan=NaN16, checkpoint
                                       )


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
savefile = get_trigger_savefile(props)
obtain_triggered_causality_binavg(em, beh, props; grd, savefile, checkpointmod=5)

hatrajcodex = Dict(r.hatrajnum=>r.hatraj for 
     r in eachrow(unique(beh[!,[:hatraj, :hatrajnum]])))

JLD2.jldsave(savefile; params, est, grd, grid_kws, 
             weighting_trig, props, ca1pfc, pfcca1, checkpoint)

# ========================
# X - Y - START - STOP
# ========================

## =====
## X - Y
## =====
#grid_kws=(;widths=[5f0,5f0])
#props = ["x","y"]
#savefile = get_trigger_savefile(props)
#ca1pfc, pfcca1, weighting_trig, checkpoint = 
#        obtain_triggered_causality_binavg(em, beh, props; grid_kws, savefile)
#
#JLD2.jldsave(savefile; params, est, grd, grid_kws, weighting_trig,
#             props, ca1pfc, pfcca1, checkpoint)

