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
using Utils.binning

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
params = (;params..., horizon=1:10, thread=false)
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

grid_kws=(;widths=[5f0,5f0])
props = ["x","y"]
propstime = ["time",props...]

EFF = EmbeddingFrameFetch(em, :area, beh, props; ordering=Dict(:area=>[:ca1,:pfc]))
grid = get_grid(beh, props; grid_kws...)
#trig = triggering.get_triggergen(0.5, grid, EFF[4]...)

ca1pfc = fill(NaN,size(grid)..., length(EFF), params[:horizon].stop)
pfcca1 = fill(NaN,size(grid)..., length(EFF), params[:horizon].stop)

#(e,eff) = first(enumerate(EFF))
@showprogress for (e,eff) in enumerate(EFF)
    zones = triggering.get_triggergen(0.5, grid, eff...)
    #(idx,(ca1,pfc)) = first(trig)
    for (idx,data) in zones
        if isempty(data); continue; end
        count = 0
        for D in data
            if isempty(D); continue; end
            ca1,pfc = map(d->transpose(hcat(d.data...)), D)
            try
            ca1pfc[idx..., e, :] .+= 
                     causal.predictiveasymmetry(ca1, pfc; params...)
            pfcca1[idx..., e, :] .+= 
                     causal.predictiveasymmetry(pfc, ca1; params...)
            catch
            end
            count += 1
        end
        ca1pfc[idx..., e, :] ./= count
        pfcca1[idx..., e, :] ./= count
    end
end

# getting trigger sample info
scatter([x for x in vec.(trig.occ.inds) if !isempty(x)]; markersize=2, markerstrokewidth=0, markerstrokealpha=0.001, markerstrokestyle=:dot, legend=:none)
