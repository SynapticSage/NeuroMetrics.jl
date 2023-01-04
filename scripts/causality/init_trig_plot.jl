#  ==================================
# PLOT TRIGGERED
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
import Munge

## ----------
## CONSTANTS
## ----------
corerr,tsk,cortsk = Munge.behavior.corerr, Munge.behavior.tsk, 
                 Munge.behavior.cortsk

## ----------
## PARAMETERS
## ----------
animal, day = "RY16", 36
esttype = :binned
est, params = get_est_preset(esttype)
params = (;params..., binning=5, window=1.25)
props = ["cuemem", "correct", "hatrajnum"]
props = ["cuemem", "correct", "hatrajnum","startWell","stopWell"]
N = 100

# Obtain savefile
# ---------------
paramstr = Utils.namedtup.tostring(params)
tagstr = if "tag" in propertynames(Main)
    "_$tag"
else
    "$animal$day.$(N)seg"
end
get_savefile = Munge.manifold.get_trigger_savefile

savefile = get_savefile(props; params)
@assert isfile(savefile)
storage = jldopen(savefile,"r");

@load savefile pfcca1
@load savefile ca1pfc
