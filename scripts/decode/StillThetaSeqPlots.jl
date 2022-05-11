# TODO: Add points ahead of and behind animal
using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
include(scriptsdir("decode","Initialize.jl"))

# Load data
include(scriptsdir("decode","LoadData.jl"))

utils.pushover("Loaded and preprocessed thetaseqplots.jl...initializing theta sequences")

@showprogress for (cyc, cycle) in collect(enumerate(eachrow(cycles)))[46237:end] # TODO hardcoding == bad!

    start, stop = cycle.start, cycle.end
    filt = (D["time"].>=start) .& (D["time"].<=stop)
    if any(filt)
        dat_sub = dat[:,:, filt]
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
    ax = Axis(fig[4,1])
    α = 0.1
    for (i,unit) in enumerate(sp)
        #print(i, " ")
        cmap = get(ColorSchemes.hawaii, ((unit.time.-start)./(stop-start)))
        cmap = parse.(Colorant, cmap)
        jitter= (α * rand(Float64,size(unit.time))) .- α/2
        try
            scatter!(ax, unit.time, i*ones(size(unit.time))*yfactor, color=cmap, strokewidth=1, markersize=6)
        catch
            @warn "$i failed"
        end
    end
    lines!(ax, lfp_sub.time, lfp_sub.raw, linewidth=2)

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
    P(t) = Point2f.(behx(t), behy(t))
    #area = cycle.area
    area = "CA1"
    color         = decode.scatter.marker_color([area])
    sc_glow_color = decode.scatter.glow_color([area])
    sc_glow_width = cycle.amp_mean/20

    # ----------------------------
    # SETUP FIGURE, AXIS, ARTISTIS
    # ----------------------------
    ax = Axis(fig[2:3,1], xlabel="centimeter", ylabel="centimeter", title="area=$(area), amp=$(round(cycle.amp_mean,digits=2))", aspect=1.7)
    boundary = Float64.(isnan.(dat_sub[:,:,1]))
    bound_color = [RGBA(colorant"grey20", 0.0), RGBA(colorant"grey14", 1.0)]
    boundary_cmap = ColorSchemeTools.make_colorscheme(bound_color, 20)
    heatmap!(ax, x, y, boundary, colormap=boundary_cmap, depth_shift=0, nan_color=RGBA(0,0,0,0))
    timestamps = range(1, size(dat_sub,3), step=1)
    cmaps = decode.heatmap.static_colormap_per_sample(:hawaii, timestamps)
    thresh = thresh_var[variable]

    for t in timestamps
        ds = decode.quantile_threshold(copy(dat_sub), thresh=thresh)
        DS = ds[:,:,t]
        hm_=heatmap!(ax, x, y, DS, overdraw=false, transparency=true, depth_shift=0+(t*0.00001),
                     colormap=cgrad(cmaps[t], alpha=0.1))
    end

    sc_ = scatter!(ax, P(1), color=color, depth_shift=1, overdraw=false,
                   transparency=true,
                   markersize=20, glowwidth=sc_glow_width,
                   glowcolor=(sc_glow_color, 0.8), label="position")
    sc_ = scatter!(Point2f.(wells.x,wells.y), color=:white,
                   strokecolor=:black, strokewidth=1.1, depth_shift=0.3,
                   overdraw=false, marker = :star5, markersize=15,
                   label="wells")

    Colorbar(fig[1, 1], limits = (0, stop-start), colormap=:hawaii,
             label = "Time", flipaxis = false, vertical=false)

    if savestuff
        for e in ["pdf","png"]
            savefile = plotsdir("theta","mpp_decode", "withBehVideo=$usevideo",
                                 outputVideo,
                                 "cycle=$cyc.$area.$tetrode.amp=$(round(cycle.amp_mean,digits=2))" *
                                 ".$e")
            save(savefile, fig, pt_per_unit = 1)
        end
    end
end

