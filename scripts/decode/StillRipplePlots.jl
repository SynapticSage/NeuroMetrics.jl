using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
include(scriptsdir("decode","Initialize.jl"))

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
beh = annotate_pastFutureGoals(beh; doPrevPast=false)

validripples = ripples[ripples.epoch.==epoch,:]
#(rip,ripple) = collect(enumerate(eachrow(validripples)))[4]
(rip,ripple) = collect(enumerate(eachrow(validripples)))[11]

@showprogress for (rip,ripple) in enumerate(eachrow(validripples))
    print(rip)

    start, stop = ripple.start, ripple.stop
    dat_sub = dat[:,:,(D["time"].>=start) .& (D["time"].<=stop)]

    if size(dat_sub,3) == 0
        continue
    end

    fig=Figure()
    B = beh[(beh.time.>=start) .& (beh.time.<=stop),:]
    sp = copy(spikes[(spikes.time.>start) .& (spikes.time .<stop),:])
    sp = groupby(sp, :unit)
    sp = [sp...]
    sp = sort(sp, by=x->median(x.time))
    axArena = Axis(fig[4,1])
    α = 0.1
    for (i,unit) in enumerate(sp)
        try
            cmap = get(ColorSchemes.hawaii, ((unit.time.-start)./(stop-start)))
            cmap = parse.(Colorant, cmap)
            jitter= (α * rand(Float64,size(unit.time))) .- α/2
            scatter!(axArena, unit.time, i*ones(size(unit.time))*4, color=cmap, strokewidth=1, markersize=3)
        catch
        end
        #if i == 1
        #    Plots.scatter(unit.time, i*ones(size(unit.time))*4)
        #else
        #    Plots.scatter!(unit.time, i*ones(size(unit.time))*4)
        #end
    end

    # -----------------
    # READY HEATMAP VARS
    # -----------------
    # -----------------------------------------------
    # Setup observable and connect decoder to behavior
    # -----------------------------------------------
    t = Int(1)
    decode_inds = SearchSortedNearest.searchsortednearest.([B.time], 
                                                           D["time"])
    # -----------------
    # READ SCATTER VARS
    # -----------------
    behx(t) = B[decode_inds[t],"x"]
    behy(t) = B[decode_inds[t],"y"]
    behpoint(t) = (behx(t), behy(t))
    now = Int(round(median(decode_inds)))
    P(t) = Point2f.(behx(t), behy(t))
    area = ripple.area
    color         = decode.scatter.marker_color([ripple.area])
    sc_glow_color = decode.scatter.glow_color([ripple.area])
    sc_glow_width = ripple.amp

    # ----------------------------
    # SETUP FIGURE, AXIS, ARTISTIS
    # ----------------------------
    axArena = Axis(fig[2:3,1], xlabel="centimeter", ylabel="centimeter", title="area=$(ripple.area), amp=$(round(ripple.amp,digits=2))")
    boundary = Float64.(isnan.(dat_sub[:,:,1]))
    bound_color = [RGBA(colorant"grey20", 0.0), RGBA(colorant"grey14", 1.0)]
    boundary_cmap = ColorSchemeTools.make_colorscheme(bound_color, 20)
    heatmap!(axArena, x, y, boundary, colormap=boundary_cmap, depth_shift=0, nan_color=RGBA(0,0,0,0))
    timestamps = range(1, size(dat_sub,3), step=1)
    cmaps = decode.heatmap.static_colormap_per_sample(:hawaii, timestamps)
    thresh = thresh_var[variable]
    for t in timestamps
        ds = decode.quantile_threshold(copy(dat_sub), thresh=thresh)
        DS = ds[:,:,t]
        hm_=heatmap!(axArena, x, y, DS, overdraw=false, transparency=true, depth_shift=0+(t*0.00001),
                     colormap=cgrad(cmaps[t], alpha=0.1))
    end
    sc_ = scatter!(axArena, P(1), color=color, depth_shift=1, overdraw=false,
                   transparency=true,
                   markersize=20, glowwidth=sc_glow_width,
                   glowcolor=(sc_glow_color, 0.8), label="position")

    sc_arena = scatter!(axArena, arenaWells.x, arenaWells.y, marker=:star5, 
                        markersize=40, color=:gray, label="Arena Well")
    sc_home  = scatter!(axArena, [homeWell.x], [homeWell.y], marker='𝐇',    
                        markersize=25, color=:gray, label="Home Well")

    future₁ = @lift(B.stopWell[now])
    future₂ = @lift(B.futureStopWell[now])
    past₁   = @lift(B.pastStopWell[now])
    correctColor = @lift begin
        if $correct == 0
            :indianred1
        elseif $correct == 1
            :mediumspringgreen
        else
            :white
        end
    end
    future₁_well_xy = @lift $future₁ == -1 ? Point2f(NaN, NaN) : Point2f(wells.x[$future₁], wells.y[$future₁])
    future₂_well_xy = @lift (($future₂ == -1) || ($future₂ == $future₁)) ? Point2f(NaN, NaN) : Point2f(wells.x[$future₂], wells.y[$future₂])
    past₁_well_xy    = @lift (($past₁ == -1) || ($past₁ == $future₁)) ? Point2f(NaN, NaN) : Point2f(wells.x[$past₁], wells.y[$past₁])
    sc_future₁_well = scatter!(axArena, future₁_well_xy,  alpha=0.5, marker='𝐅', markersize=60, color=correctColor, glowwidth=29)
    sc_future₂_well = scatter!(axArena, future₂_well_xy,  alpha=0.5, marker='𝐟', markersize=60, color=:white,       glowwidth=5)
    sc_past_well    = scatter!(axArena, past₁_well_xy,    alpha=0.5, marker='𝐏', markersize=60, color=:white,       glowwidth=5)
    if doPrevPast
        past₂   = @lift($behavior.pastStopWell[now])
        past₂_well_xy    = @lift (($past₂ == -1) || ($past₂ == $future₁)) ? Point2f(NaN, NaN) : Point2f(wells.x[$past₂], wells.y[$past₂])
    end

    Colorbar(fig[1, 1], limits = (0, stop-start), colormap=:hawaii,
             label = "Time", flipaxis = false, vertical=false)
    if savestuff
        for e in ["pdf","svg"]
            savefile = plotsdir("ripples","mpp_decode", "withBehVideo=$usevideo",
                                 outputVideo,
                                 "rip=$rip.$area.amp=$(round(ripple.amp,digits=2))" *
                                 ".$e")
            save(savefile, fig, pt_per_unit = 1)
        end
    end

end

