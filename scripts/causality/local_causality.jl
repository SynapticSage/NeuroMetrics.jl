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
lab = OrderedDict([0,1]=>"CUE correct", [1,1]=>"MEM correct",
                  [0,0]=>"CUE error", [1,0]=>"MEM error")
conditionals = [:cuemem,:correct]
animal, day = "RY16", 36

# Obtain params
# ---------------
@info "prop missinG"
esttype = :binned
est, params = get_est_preset(esttype)
params = (;params..., horizon=1:20, thread=false)
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

savefile = datadir("manifold","causal","local_grid_cause_props=$(join(props,","))_$(paramstr)_$tagstr.jld2")
if isfile(savefile)
    D=JLD2.load(savefile)
    Utils.dict.load_dict_to_module!(Main, D)
end

@assert !isempty(embedding)
area_embeddings = groupby(em, :area)

grid_kws=(;widths=[5f0,5f0])
props = ["x","y"]
propstime = ["time",props...]

EFF = EmbeddingFrameFetch(em, :area, beh, props; 
                          ordering=Dict(:area=>[:ca1,:pfc]))
grid = get_grid(beh, props; grid_kws...)
#trig = triggering.get_triggergen(0.5, grid, EFF[4]...)

ca1pfc = fill(NaN,size(grid)..., length(EFF), params[:horizon].stop)
pfcca1 = fill(NaN,size(grid)..., length(EFF), params[:horizon].stop)

#(e,eff) = first(enumerate(EFF))
@showprogress "frames" for (e,eff) in enumerate(EFF[2:length(EFF)])
    @info e
    zones = triggering.get_triggergen(1, grid, eff...)
    #(idx,(ca1,pfc)) = first(trig)
    @showprogress "zones" for (idx,data) in zones
        if isempty(data); continue; end
        count = 0
        cpv = view(ca1pfc, idx..., e, :)
        pcv = view(pfcca1, idx..., e, :)
        for D in data
            if isempty(D); continue; end
            ca1,pfc = map(d->transpose(hcat(d.data...)), D)
            try
                cpt = causal.predictiveasymmetry(ca1, pfc; params...)
                pct = causal.predictiveasymmetry(pfc, ca1; params...)
                if all(isnan.(cpv)) && all(isnan.(pcv))
                    cpv .= 0
                    pcv .= 0
                end
                cpv .+= cpt
                pcv .+= pct
            catch err
                showerror(stdout, err)
                @info "sizes" size(ca1) size(pfc)
                @infiltrate
            end
            count += 1
        end
        @infiltrate
        ca1pfc[idx..., e, :] ./= count
        pfcca1[idx..., e, :] ./= count
    end
end

JLD2.jldsave(savefile; params, est, grid, grid_kws, 
             props, propstime, ca1pfc, pfcca1)

# How many unsampled?
heatmap((sum(isnan.(ca1pfc), dims=(3,4))./size(ca1pfc,4))[:,:])

# getting trigger sample info
#scatter([x for x in vec.(trig.occ.inds) if !isempty(x)]; markersize=2, markerstrokewidth=0, markerstrokealpha=0.001, markerstrokestyle=:dot, legend=:none)
C = [collect(Float64.(c)) for c in grid.centers]
h1=heatmap(C[1], C[2], nanmean(pfcca1, dims=(3,4))[:,:]', c=:vik, clims=(-0.006,0.006))
plotboundary!(tsk,transpose=true)
h2=heatmap(C[1], C[2], nanmean(ca1pfc, dims=(3,4))[:,:]', c=:vik, clims=(-0.006,0.006))
plotboundary!(tsk,transpose=true)
plot(h1,h2, size=(800,400))

