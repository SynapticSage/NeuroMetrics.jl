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


## ----------
## PARAMETERS
## ----------

function obtain_local_binned_measure(em, beh, props; grid_kws=nothing, savefile,
        grd=(println("default binning compute");binning.get_grid(beh, props; grid_kws...)),
        checkpoint=ThreadSafeDict{NamedTuple, Bool}())
    if isfile(savefile)
        D=JLD2.load(savefile)
        Utils.dict.load_dict_to_module!(Main, D)
    end
    ## ------------------------------
    ## Setup frames and zone triggers
    ## ------------------------------
    propstime = ["time",props...]
    EFF = EmbeddingFrameFetch(em, :area, beh, props; 
                              ordering=Dict(:area=>[:ca1,:pfc]))

    ## ------------------------
    ## Setup stores for results
    ## ------------------------
    ca1pfc = fill(NaN,size(grd)..., length(EFF), params[:horizon].stop)
    pfcca1 = fill(NaN,size(grd)..., length(EFF), params[:horizon].stop)
    

    ## ------------------------
    ## Compute
    ## ------------------------
    #(e,eff) = first(enumerate(EFF))

    @info "How much ground to cover" length(EFF) length(first(EFF))

    @showprogress "frames" for (e,eff) in enumerate(EFF[1:length(EFF)])
        @info e
        zones = triggering.get_triggergen(1, grd, eff...)
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
            for D in data
                if isempty(D); continue; end
                ca1, pfc = map(d->transpose(hcat(d.data...)), D)
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
                    #@infiltrate
                end
                count += 1
            end
            ca1pfc[idx..., e, :] ./= count
            pfcca1[idx..., e, :] ./= count
            checkpoint[checkpointkey] = true
        end
        storage    = JLD2.jldopen(savefile, "a")
        delete!(storage,"ca1pfc")
        storage["ca1pfc"] = ca1pfc
        delete!(storage,"pfcca1")
        storage["pfcca1"] = pfcca1
        delete!(storage,"checkpoint")
        storage["checkpoint"] = checkpoint
        close(storage)
    end


    # Better weighting
    weighting_trig = zeros(size(grd)..., length(EFF))
    times = []
    for (e,eff) in enumerate(EFF)
        zones = triggering.get_triggergen(1, grd, eff...)
        #@infiltrate length(zones) == 0
        for (idx,data) in zones
            weighting_trig[idx..., e] += length(data)
            [append!(times, d[1].time) for d in data]
            #@infiltrate !isempty(data)
        end
        # @infiltrate all(weighting_trig[:,:,e] .== 0)
    end
    (;ca1pfc, pfcca1, weighting_trig, checkpoint)
end
get_savefile(props) = datadir("manifold","causal","local_grid_cause_props=$(join(props,","))_$(paramstr)_$tagstr.jld2")

# =====
# X - Y
# =====
grid_kws=(;widths=[5f0,5f0])
props = ["x","y"]
savefile = get_savefile(props)
ca1pfc, pfcca1, weighting_trig, checkpoint = 
        obtain_local_binned_measure(em, beh, props; grid_kws, savefile)
JLD2.jldsave(savefile; params, est, grd, grid_kws, weighting_trig,
             props, ca1pfc, pfcca1, checkpoint)

# ========================
# X - Y - CUEMEM - CORRERR
# ========================
props = ["x", "y", "cuemem", "correct"]
grid_kws=(;widths=[5f0,5f0,1f0,1f0], radiusinc=[1.2f0,1.2f0,1f0,1f0])
grd = Utils.binning.get_grid(beh, props; grid_kws...)
savefile = get_savefile(props)
ca1pfc, pfcca1, weighting_trig, checkpoint = 
        obtain_local_binned_measure(em, beh, props; grd, savefile)

# ========================
# X - Y - START - STOP
# ========================

#                    
#,---.|         |    
#|---'|    ,---.|--- 
#|    |    |   ||    
#`    `---'`---'`---'
#                    

func = nanmean

using Plot.task
tsk = Load.load_task(animal, day)

C = [collect(Float64.(c)) for c in grd.centers]
# How many unsampled?
nanfrac = (sum(isnan.(ca1pfc), dims=(3,4))./size(ca1pfc,4))[:,:]'./size(ca1pfc,3)
S=heatmap(C..., nanfrac, title="unsampled fraction")
plotboundary!(tsk,transpose=true, c=:black)

weightstyle=:nanfrac
if weightstyle == :nanfrac
    weighting = (1 .- nanfrac)'
elseif weightstyle == :sampcount
    weighting = weighting_trig
elseif weightstyle == :nanfrac_sq
    weighting = ((1 .- nanfrac) .^ 2)'
end

function plot_track_cause(ind=Colon())
    # getting trigger sample info
    #scatter([x for x in vec.(trig.occ.inds) if !isempty(x)]; markersize=2, markerstrokewidth=0, markerstrokealpha=0.001, markerstrokestyle=:dot, legend=:none)
    #heatmap(nanmean(pfcca1 - ca1pfc, dims=(3,4))[:,:])
    ind = ind isa Real ? [ind] : ind
    X=func(pfcca1[:,:,:,ind] .* weighting, dims=(3,4))[:,:]' 
    h1=heatmap(C[1], C[2], X, title="pfc-ca1", c=:vik, clims=nanmaximum(abs.(X)).*(-1,1))
    plotboundary!(tsk,transpose=true, c=:black)
    X=func(ca1pfc[:,:,:,ind] .* weighting, dims=(3,4))[:,:]' 
    h2=heatmap(C[1], C[2], X, title="ca1-pfc", c=:vik, clims=nanmaximum(abs.(X)).*(-1,1))
    plotboundary!(tsk,transpose=true, c=:black)
    plot(h1,h2, size=(800,400))

    plot(h1,h2,S, size=(800,800))

    dif = pfcca1 .- ca1pfc
    dif = func(dif[:,:,:,ind] .* weighting, dims=(3,4))[:,:]' .* (2)
    hD = heatmap(C..., dif; title="pfcca1 .- ca1pfc", c=:vik, clims=nanmaximum(abs.(dif)).*(-1,1))
    plotboundary!(tsk,transpose=true,c=:black,linestyle=:solid,label="")

    plot(h1,h2,S,hD, size=(800,400))
end
@time plot_track_cause()

anim = @animate for i in params[:horizon]
    plot_track_cause(i)
end
gif(anim, fps=5)
