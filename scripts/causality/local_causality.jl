using DrWatson
using GoalFetchAnalysis
using Infiltrator
import Plot
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

cor = Dict(0=>"correct", 1=>"error")
tsk = Dict(0=>"cue", 1=>"mem")
lab = OrderedDict([0,1]=>"CUE correct", [1,1]=>"MEM correct", [0,0]=>"CUE error", [1,0]=>"MEM error")
conditionals = [:cuemem,:correct]
animal, day = "RY16", 36

# Obtain params
# ---------------
@info "prop missinG"
esttype = :binned
est, params = get_est_preset(esttype)
params = (;params..., horizon=1:250, thread=true)
params = (;params..., binning=5)

# Obtain savefile
# ---------------
paramstr = Utils.namedtup.tostring(params)
N = 100
tagstr = if "tag" in propertynames(Main)
    "_$tag"
else
    "$animal$day.$(N)seg"
end
savefile = datadir("manifold","causal","pa_cause_$(paramstr)_$tagstr.jld2")


D=JLD2.load(savefile)
Utils.dict.load_dict_to_module!(Main, D)

@assert !isempty(embedding)
area_embeddings = groupby(em, :area)

EFF = EmbeddingFrameFetch(CA1PFC, beh, ["x","y"])
trig = triggering.get_triggergen(, ["x","y"], 0.5; widths=[5f0,5f0])

