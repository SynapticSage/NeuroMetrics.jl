# TODO: Add points ahead of and behind animal
using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
include(scriptsdir("decode","InitializeCairo.jl"))

for (splitfig, split_num) in Iterators.product([true,false], 0:3)

    thresh = Dict("likelihood"=>0.1, "acausal_posterior"=>0.9875, 
                  "causal_posterior"=> 0.9875)
    dothresh  = false
    dodisplay = true
    spikecolorby = :time # :celltime | :time
    padding = [.10, .10]
    @info split_num
    global decode_file
    global beh, spikes, cycles, dat, T, x, y
    function get_color(x,cmap,nan_color,alpha_all=1)
        if isnan(x)
            nan_color
        else
            col=get(cgrad(cmap), x)
            RGBA(col.r, col.g, col.b, alpha_all)
        end
    end


    # Load data
    #include(scriptsdir("decode","Initialize.jl"))
    include(scriptsdir("decode","InitializeCairo.jl"))
    decode_file=replace(decode_file, "split=0"=>"split=$split_num")
    @info decode_file
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
    (cyc ,cycle) = collect(enumerate(eachrow(cycles)))[2]

    @showprogress for cycle in collect(eachrow(cycles)) # TODO hardcoding == bad!

        cyc = cycle.cycle
        @info cyc

        start, stop = cycle.start, cycle.stop
        filt_cyc = (T.>=start) .& (T.<=stop)
        dat_sub = dat[:,:, filt_cyc]
        lfp_filt = (lfp[!,"time"].>=(start-padding[1])) .&
                   (lfp[!,"time"].<=(stop +padding[2]))
        lfp_sub = lfp[lfp_filt, :]
        lfp_sub_filt = (lfp_sub[!,"time"].>=start) .&
                       (lfp_sub[!,"time"].<=stop)

        if !(any(filt_cyc)) || isempty(dat_sub)
            continue
        end
    sp_filt = (spikes.time .>= (start-padding[1])) .&
              (spikes.time .<  (stop+padding[2]))
    sp = copy(spikes[sp_filt,:])
    sp.reltime = max.(min.((sp.time .- start)./(stop-start), 1), 0)
    sp.ranreltime = (sp.time .- start)./(stop-start)
    #TODO splpit by cell andd do this onlly when no spikes in ripple
    sp.ranreltime[sp.ranreltime .< 0 .||  sp.ranreltime .> 1] = Random.shuffle(sp.ranreltime[sp.ranreltime .< 0 .||  sp.ranreltime .> 1])

    future₁ = B.stopWell[1]
    future₂ = B.futureStopWell[1]
    past₁   = B.pastStopWell[1]


        B = beh[(beh.time.>=start) .& (beh.time.<=stop),:]
        if isempty(B)
            @warn "Empty behavior for cycle $cyc"
            continue
        end
        area = "CA1"

        if !(splitfig)
            fig = Figure(resolution=Tuple(0.8*[1200,1600]))
            axArena = Axis(fig[2:3,1], xlabel="centimeter", ylabel="centimeter",
                           title="area=$(area), amp=$(round(cycle.amp_mean,digits=2))\nfuture₁=$future₁, future₂=$future₂,past₁=$past₁", aspect=1.7)
            axSpikes = Axis(fig[4,1])
            Colorbar(fig[1, 1], limits = (0, stop-start), colormap=:hawaii,
                     label = "Time", flipaxis = false, vertical=false)
            bc = fig.scene.backgroundcolor[]
        else
            figSpikes = Figure(resolution=Tuple(0.8*[1200,800]))
            figArena  = Figure(resolution=Tuple(0.8*[1200,800]))
            axArena = Axis(figArena[2:3,1], xlabel="centimeter", ylabel="centimeter",
                           title="area=$(area), amp=$(round(cycle.amp_mean,digits=2))\nfuture₁=$future₁, future₂=$future₂,past₁=$past₁", aspect=1.7)
            axSpikes = Axis(figSpikes[1,1])
            Colorbar(figArena[1, 1], limits = (0, stop-start), colormap=:hawaii,
                     label = "Time", flipaxis = false, vertical=false)
            bc = figSpikes.scene.backgroundcolor[]
        end

        yfactor=4
        lfp_sub.raw = utils.norm_extrema(lfp_sub.raw, [0, length(unique(sp.unit))*yfactor])
        sp = groupby(sp, :unit)
        sp = [sp...]
        sp = sort(sp, by=x->median(x.reltime))
        α = 0.1
        for (i,unit) in enumerate(sp)
            inside_range = .&(unit.time .>= start,
                              unit.time .<  stop)
            if spikecolorby == :celltime
                med_utime = (median(unit.time).-start)./(stop-start)
                reltime = max.(min.(med_utime, 1),0)
                color = get(ColorSchemes.hawaii, reltime)
                color = repeat([color], length(unit.time))
            elseif spikecolorby == :time
                med_utime = (unit.time.-start)./(stop-start)
                reltime   = max.(min.(med_utime, 1),0)
                color = get(ColorSchemes.hawaii, reltime)
            end
            nan_color = RGBA(1,1,1,1)
            color[(!).(inside_range)] .= nan_color
            try
                scatter!(axSpikes, unit.time, i*ones(size(unit.time))*yfactor,
                         color=color, strokewidth=1, markersize=7)
            catch
                @warn "$i failed"
            end
        end
        nan_color = RGBA(1,1,1,1)
        zerod_phase = ((lfp_sub.time .-start))./(stop-start)
        zerod_phase[(zerod_phase .> 1) .| (zerod_phase .< 0)] .= NaN
        zerod_phase = get_color.(zerod_phase, :hawaii, nan_color)
        lines!(axSpikes, lfp_sub.time, lfp_sub.raw, colormap=:hawaii,
               color=zerod_phase, linewidth=2)

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
        color         = decode.scatter.marker_color([area])
        sc_glow_color = decode.scatter.glow_color([area])
        sc_glow_width = cycle.amp_mean/20

        # ----------------------------
        # SETUP FIGURE, AXIS, ARTISTIS
        # ----------------------------
        lines!(axArena, boundary.x, boundary.y, color=:grey)
        timestamps = range(1, size(dat_sub,3), step=1)
        cmaps = decode.heatmap.static_colormap_per_sample(:hawaii, timestamps)
        cmaps = Vector{Any}([cmaps...])
        thr= thresh[variable]
        ds = decode.quantile_threshold(copy(dat_sub), thr)

        backgroundcolor = RGBA(bc.r, bc.g, bc.b, bc.alpha)
        nan_color = RGBA(bc.r, bc.g, bc.b, 0)
        DS = cat([get_color.(ds[:,:,t], [cmaps[t]], nan_color, 0.45) for t in 1:size(ds,3)]...;
                 dims=3)

        [heatmap!(axArena, x, y, DS[:,:,t], transparency=true, overdraw=false,
                  depth_shift=0+(t*0.00001)) for t in 1:size(DS,3)]

        sc_ = scatter!(axArena, P(1), color=color, depth_shift=1, overdraw=false,
                       transparency=true,
                       markersize=20, glowwidth=sc_glow_width,
                       glowcolor=(sc_glow_color, 0.8), label="position")
        sc_arena = scatter!(axArena, arenaWells.x, arenaWells.y, marker=:star5, 
                            markersize=40, color=:gray, label="Arena Well")
        sc_home  = scatter!(axArena, [homeWell.x], [homeWell.y], marker='𝐇',    
                            markersize=25, color=:gray, label="Home Well")
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
        past₁_well_xy   = ((past₁ == -1) || (past₁ == future₁)) ? Point2f(NaN, NaN) : Point2f(wells.x[past₁], wells.y[past₁])
        sc_future₁_well = scatter!(axArena, future₁_well_xy,  alpha=0.5, marker='𝐅', markersize=60, color=correctColor, glowwidth=29)
        sc_future₂_well = scatter!(axArena, future₂_well_xy,  alpha=0.5, marker='𝐟', markersize=60, color=:white,       glowwidth=5)
        sc_past_well    = scatter!(axArena, past₁_well_xy,    alpha=0.5, marker='𝐏', markersize=60, color=:white,       glowwidth=5)
        if doPrevPast
            past₂         = B.pastStopWell[1]
            past₂_well_xy = ((past₂ == -1) || (past₂ == future₁)) ? Point2f(NaN, NaN) : Point2f(wells.x[past₂], wells.y[past₂])
        end

        if dodisplay && !(splitfig)
           electrondisplay(fig)
        elseif dodisplay && splitfig
           electrondisplay(figSpikes)
           electrondisplay(figArena)
        end

        if savestuff
            for e in ["pdf"]
                if splitfig
                    spikesFile = plotsdir("theta","mpp_decode", "withBehVideo=$usevideo", outputVideo,
                                         "spikes=$cyc.$area.$tetrode.amp=$(round(cycle.amp_mean,digits=2))" *
                                         ".$e")
                    save(spikesFile, figSpikes, pt_per_unit=1)
                    arenaFile = plotsdir("theta","mpp_decode", "withBehVideo=$usevideo", outputVideo,
                                         "arena=$cyc.$area.$tetrode.amp=$(round(cycle.amp_mean,digits=2))" *
                                         ".$e")
                    save(arenaFile, figArena, pt_per_unit=1)
                else
                    savefile = plotsdir("theta","mpp_decode", "withBehVideo=$usevideo", outputVideo,
                                         "cycle=$cyc.$area.$tetrode.amp=$(round(cycle.amp_mean,digits=2))" * ".$e")
                    save(savefile, fig, pt_per_unit = 1)
                end
            end
        end

    end
end
