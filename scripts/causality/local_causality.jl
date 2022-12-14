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

function obtain_local_binned_measure(em, beh, props; grid_kws=nothing, savefile,
        grd=(println("default binning compute");binning.get_grid(beh, props; grid_kws...)),
        checkpoint=ThreadSafeDict{NamedTuple, Bool}())
    if isfile(savefile)
        D=JLD2.load(savefile)
        Utils.dict.load_dict_to_module!(Main, D)
    end
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
    ca1pfc = fill(NaN,size(grd)..., length(EFF), params[:horizon].stop)
    pfcca1 = fill(NaN,size(grd)..., length(EFF), params[:horizon].stop)
    

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
        for (idx,data) in zones
            weighting_trig[idx..., e] += length(data)
            [append!(times, d[1].time) for d in data]
            #@infiltrate !isempty(data)
        end
        # @infiltrate all(weighting_trig[:,:,e] .== 0)
    end
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
        obtain_local_binned_measure(em, beh, props; grid_kws, savefile)

JLD2.jldsave(savefile; params, est, grd, grid_kws, weighting_trig,
             props, ca1pfc, pfcca1, checkpoint)
# ==========================
# CUEMEM - CORRECT - HA_TRAJ
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
        obtain_local_binned_measure(em, beh, props; grd, savefile)

JLD2.jldsave(savefile; params, est, grd, grid_kws, weighting_trig,
             props, ca1pfc, pfcca1, checkpoint)


# ========================
# X - Y - CUEMEM - CORRERR
# ========================
props = ["x", "y", "cuemem", "correct", "hatrajnum"]
grid_kws =
        (;widths    = [4f0,4f0,1f0,1f0,1f0],
          radiusinc = [0.2f0,0.2f0,0f0,0f0,0f0],
          maxrad    = [6f0,6f0,0.5f0,0.5f0],
         )
grd = Utils.binning.get_grid(beh, props; grid_kws...)
savefile = get_savefile(props)
ca1pfc, pfcca1, weighting_trig, checkpoint = 
        obtain_local_binned_measure(em, beh, props; grd, savefile)

JLD2.jldsave(savefile; params, est, grd, grid_kws, weighting_trig,
             props, ca1pfc, pfcca1, checkpoint)

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

function getnanfrac(X)
    horizonAndSamp = Tuple(collect(1:ndims(X))[end-1:end])
    total_who_could_be_present = size(X,horizonAndSamp[1]) .* 
                                 size(X,horizonAndSamp[2])
    sum(isnan.(X); dims=horizonAndSamp) ./
                        total_who_could_be_present
end
"""
    plottrackcause(ca1pfc,pfcca1, 
                     slice=[Colon() for _ in 1:ndims(ca1pfc)], 
                          plotdims=1:2; weightstyle=:nanfrac)

Plot the causal data in ca1pfc and pfcca1 matrices, indexed by `slice`, and
along the plotdimensions `plotdims`
"""
function plottrackcause(ca1pfc,pfcca1, 
                          slice=[Colon() for _ in 1:ndims(ca1pfc)], 
                          plotdims=1:2; weightstyle=:nanfrac)

    axes = [collect(Float64.(c)) for c in grd.centers][plotdims]

    # How many unsampled?
    nanfrac = getnanfrac(ca1pfc)

    # Perform slicing
    ca1pfc = ca1pfc[slice...]
    pfcca1 = pfcca1[slice...]
    nanfrac = nanfrac[slice...]

    squishdims = (setdiff(1:ndims(ca1pfc), collect(plotdims))...,)
    ca1pfc  = nanmean(ca1pfc,dims=squishdims)
    pfcca1  = nanmean(pfcca1;dims=squishdims)
    nanfrac = nansum(nanfrac;dims=squishdims)

    # Make sure the final product is XD
    ca1pfc  = Utils.squeeze(ca1pfc)
    pfcca1  = Utils.squeeze(pfcca1)
    nanfrac = Utils.squeeze(nanfrac)

    if weightstyle == :nanfrac
        weighting = (1 .- nanfrac)
    elseif weightstyle == :sampcount
        weighting = weighting_trig'
    elseif weightstyle == :nanfrac_sq
        weighting = ((1 .- nanfrac) .^ 2)
    end
    weigthing = replace(weighting, 0.0 => NaN)

    pfcca1  = (pfcca1 .* weighting)[:,:]'
    ca1pfc  = (ca1pfc .* weighting)[:,:]'
    nanfrac = nanfrac'

    @assert ndims(ca1pfc)<=2 "Dafuq. Dims need to be no more than 2."
    
    p_unsamp=heatmap(axes..., nanfrac, title="unsampled fraction")
    plotboundary!(tsk,transpose=true, c=:black)

    
    h_pfcca1 = heatmap(axes[1], axes[2], pfcca1, title="pfc-ca1", 
                     c=:vik, clims=nanmaximum(abs.(pfcca1)).*(-1,1))
    plotboundary!(tsk,transpose=true, c=:black)

    h_ca1pfc=heatmap(axes[1], axes[2], ca1pfc, title="ca1-pfc", 
                     c=:vik, clims=nanmaximum(abs.(ca1pfc)).*(-1,1))
    plotboundary!(tsk,transpose=true, c=:black)

    plot(h_pfcca1,h_ca1pfc, size=(800,400))
    plot(h_pfcca1,h_ca1pfc,p_unsamp, size=(800,800))

    dif = pfcca1 .- ca1pfc
    dif = func(dif .* weighting, dims=(3,4))[:,:]' .* (2)
    hD = heatmap(axes..., dif; title="pfcca1 .- ca1pfc", c=:vik, clims=nanmaximum(abs.(dif)).*(-1,1))
    plotboundary!(tsk,transpose=true,c=:black,linestyle=:solid,label="")

    plot(h_pfcca1,h_ca1pfc,p_unsamp,hD, size=(800,400))
end
@time plottrackcause()

anim = @animate for i in params[:horizon]
    plottrackcause(i)
end
gif(anim, fps=5)
