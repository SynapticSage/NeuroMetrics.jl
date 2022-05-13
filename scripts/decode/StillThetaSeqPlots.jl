# TODO: Add points ahead of and behind animal
using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
#include(scriptsdir("decode","Initialize.jl"))
include(scriptsdir("decode","InitializeCairo.jl"))
thresh = Dict("likelihood"=>0.1, "acausal_posterior"=>0.9875, 
              "causal_posterior"=> 0.9875)
dothresh=false
dodisplay=false

# Load data
include(scriptsdir("decode","LoadData.jl"))
@time task   = raw.load_task(animal,     day)
wells = task[(task.name.=="welllocs") .& (task.epoch .== epoch), :]
homeWell, arenaWells = begin
    hw = argmax(StatsBase.fit(StatsBase.Histogram, filter(b->b!=-1,beh.stopWell), 1:6).weights)
    @info "Homewell = $hw"
    wells[hw,:], wells[setdiff(1:5,hw),:]
end
boundary = task[(task.name.=="boundary") .& (task.epoch .== epoch), :]
append!(boundary, DataFrame(boundary[1,:]))
include(scriptsdir("decode", "PreprocessLFP.jl"))
beh = annotate_pastFutureGoals(beh; doPrevPast=false)
utils.pushover("Loaded and preprocessed thetaseqplots.jl...initializing theta sequences")
savestuff = true
tetrode   = 5

@showprogress for (cyc, cycle) in collect(enumerate(eachrow(cycles))) # TODO hardcoding == bad!

    @info "cyc=$cyc"
    start, stop = cycle.start, cycle.stop
    filt_cyc = (T.>=start) .& (T.<=stop)
    if any(filt_cyc)
        dat_sub = dat[:,:, filt_cyc]
        lfp_filt = (lfp[!,"time"].>=start) .& (lfp[!,"time"].<=stop)
        lfp_sub = lfp[lfp_filt, :]
        #println("cycle ", cyc)
        #break
    else
        continue
    end

    B = beh[(beh.time.>=start) .& (beh.time.<=stop),:]
    if isempty(B)
        @warn "Empty behavior for cycle $cyc"
        continue
    end

    fig=Figure(resolution=Tuple(0.8*[1200,1600]))
    yfactor=4
    sp = copy(spikes[(spikes.time.>start) .& (spikes.time .<stop),:])
    lfp_sub.raw = utils.norm_extrema(lfp_sub.raw, [0, length(unique(sp.unit))*yfactor])
    sp = groupby(sp, :unit)
    sp = [sp...]
    sp = sort(sp, by=x->median(x.time))
    axArena = Axis(fig[4,1])
    α = 0.1
    for (i,unit) in enumerate(sp)
        #print(i, " ")
        cmap = get(ColorSchemes.hawaii, ((unit.time.-start)./(stop-start)))
        cmap = parse.(Colorant, cmap)
        try
            scatter!(axArena, unit.time, i*ones(size(unit.time))*yfactor, color=cmap, strokewidth=1, markersize=6)
        catch
            @warn "$i failed"
        end
    end
    lines!(axArena, lfp_sub.time, lfp_sub.raw, linewidth=2)

    # -----------------
    # READY HEATMAP VARS
    # -----------------
    # -----------------------------------------------
    # Setup observable and connect decoder to behavior
    # -----------------------------------------------
    t = Int(1)
    decode_inds = utils.searchsortednearest.([B.time], T)

    # -----------------
    # READ SCATTER VARS
    # -----------------
    behx(t) = B[decode_inds[t],"x"]
    behy(t) = B[decode_inds[t],"y"]
    behpoint(t) = (behx(t), behy(t))
    P(t) = Point2f.(behx(t), behy(t))
    now = Int(round(median(decode_inds)))
    #area = cycle.area
    area = "CA1"
    color         = decode.scatter.marker_color([area])
    sc_glow_color = decode.scatter.glow_color([area])
    sc_glow_width = cycle.amp_mean/20

    # ----------------------------
    # SETUP FIGURE, AXIS, ARTISTIS
    # ----------------------------
    future₁ = B.stopWell[1]
    future₂ = B.futureStopWell[1]
    past₁   = B.pastStopWell[1]
    axArena = Axis(fig[2:3,1], xlabel="centimeter", ylabel="centimeter",
                   title="area=$(area), amp=$(round(cycle.amp_mean,digits=2))\nfuture₁=$future₁, future₂=$future₂,past₁=$past₁", aspect=1.7)
    lines!(axArena, boundary.x, boundary.y, color=:grey)
    timestamps = range(1, size(dat_sub,3), step=1)
    cmaps = decode.heatmap.static_colormap_per_sample(:hawaii, timestamps)
    cmaps = Vector{Any}([cmaps...])
    thr= thresh[variable]
    ds = decode.quantile_threshold(copy(dat_sub), thr)

    bc = fig.scene.backgroundcolor[]
    backgroundcolor = RGBA(bc.r, bc.g, bc.b, bc.alpha)
    nan_color = RGBA(bc.r, bc.g, bc.b, 0)
    function get_color(x,cmap)
        if isnan(x)
            nan_color
        else
            get(cgrad(cmap), x)
        end
    end

    DS = cat([get_color.(ds[:,:,t], [cmaps[t]]) for t in 1:size(ds,3)]...;
             dims=3)

    [heatmap!(axArena, x, y, DS[:,:,t], transparency=true, overdraw=false, depth_shift=0+(t*0.00001))
     for t in 1:size(DS,3)]

    sc_ = scatter!(axArena, P(1), color=color, depth_shift=1, overdraw=false,
                   transparency=true,
                   markersize=20, glowwidth=sc_glow_width,
                   glowcolor=(sc_glow_color, 0.8), label="position")
    sc_arena = scatter!(axArena, arenaWells.x, arenaWells.y, marker=:star5, 
                        markersize=40, color=:gray, label="Arena Well")
    sc_home  = scatter!(axArena, [homeWell.x], [homeWell.y], marker='𝐇',    
                        markersize=25, color=:gray, label="Home Well")
    Colorbar(fig[1, 1], limits = (0, stop-start), colormap=:hawaii,
             label = "Time", flipaxis = false, vertical=false)
    correctColor = begin
        if B.correct[1] == 0
            :indianred1
        elseif B.correct[1] == 1
            :mediumspringgreen
        else
            :white
        end
    end
    future₁_well_xy = future₁ == -1 ? Point2f(NaN, NaN) : Point2f(wells.x[future₁], wells.y[future₁])
    future₂_well_xy = ((future₂ == -1) || (future₂ == future₁)) ? Point2f(NaN, NaN) : Point2f(wells.x[future₂], wells.y[future₂])
    past₁_well_xy    = ((past₁ == -1) || (past₁ == future₁)) ? Point2f(NaN, NaN) : Point2f(wells.x[past₁], wells.y[past₁])
    sc_future₁_well = scatter!(axArena, future₁_well_xy,  alpha=0.5, marker='𝐅', markersize=60, color=correctColor, glowwidth=29)
    sc_future₂_well = scatter!(axArena, future₂_well_xy,  alpha=0.5, marker='𝐟', markersize=60, color=:white,       glowwidth=5)
    sc_past_well    = scatter!(axArena, past₁_well_xy,    alpha=0.5, marker='𝐏', markersize=60, color=:white,       glowwidth=5)
    if doPrevPast
        past₂   = B.pastStopWell[1]
        past₂_well_xy    = ((past₂ == -1) || (past₂ == future₁)) ? Point2f(NaN, NaN) : Point2f(wells.x[past₂], wells.y[past₂])
    end

    if dodisplay
        electrondisplay(fig)
    end

    if savestuff
        for e in ["pdf",]
            savefile = plotsdir("theta","mpp_decode", "withBehVideo=$usevideo",
                                 outputVideo,
                                 "cycle=$cyc.$area.$tetrode.amp=$(round(cycle.amp_mean,digits=2))" *
                                 ".$e")
            save(savefile, fig, pt_per_unit = 1)
        end
    end

end

