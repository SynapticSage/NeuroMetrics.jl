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
corerr,tsk,lab = Munge.behavior.corerr, Munge.behavior.tsk, 
                Munge.behavior.cortsk

## ----------
## PARAMETERS
## ----------
animal, day = "RY16", 36
esttype = :binned
est, params = get_est_preset(esttype)
params = (;params..., binning=7, window=1.25)
params = (;params..., horizon=1:30, thread=false)
props = ["cuemem", "correct", "hatrajnum"]
props = ["cuemem", "correct", "hatrajnum","startWell","stopWell"]
N = 100
spikes, beh, ripples, cells  = Load.load(animal, day)

# ---------------
# Obtain savefile
# ---------------
storage = load_alltimes_savefile(animal, day, N; params)
