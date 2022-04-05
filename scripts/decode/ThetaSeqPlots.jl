using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
using Pushover, Revise, Interact, Blink, Mux, ProgressMeter
using Statistics
using VideoIO
using ColorSchemes, Colors
import ColorSchemeTools
import DSP
savestuff = true
if savestuff
    using CairoMakie, Makie
else
    using GLMakie, Makie
end
using DataFrames
import Plots

set_theme!(theme_dark())
__revise_mode__ = :eval
#includet(srcdir("table.jl"))
includet(srcdir("raw.jl"))
includet(srcdir("utils.jl"))
includet(srcdir("table.jl"))
includet(srcdir("raster.jl"))
includet(srcdir("decode.jl"))
includet("/home/ryoung/Code/projects/goal-code/src/utils/SearchSortedNearest.jl/src/SearchSortedNearest.jl")

# -----------
# DATASETS
# -----------
animal, day, epoch = "RY16", 36, 7
@time beh    = raw.load_behavior(animal, day)
@time spikes = raw.load_spikes(animal,   day)
@time cells = raw.load_cells(animal,   day)
@time tetrode = raw.load_tetrode(animal,   day)
@time lfp  = raw.load_lfp(animal,      day)
@time task   = raw.load_task(animal,     day)
utils.pushover("Finished loaded thetaseqplots.jl")

# -----------
# SETTINGS
# -----------
usevideo=false
transition_type ="empirical"
decoder_type ="sortedspike"
variable = "causal_posterior"
video="/Volumes/Colliculus/RY16_experiment/actualVideos/$(animal)_69_0$(epoch)_CMt.1.mp4"
(split, split_type) = collect(Iterators.product([0,1,2], ["test","train"]))[1]
thresh_var = Dict("likelihood"=>0.1,
                  "acausal_posterior"=>0.985,
                  "causal_posterior"=>0.985)
outputVideo = "animation.$(decoder_type)_$(transition_type)_$(split_type)_$(split)_$(variable)_$(basename(video))"
utils.mkifne(plotsdir("ripples","mpp_decode", "withBehVideo=$usevideo", outputVideo))
utils.mkifne(plotsdir("ripples","mpp_decode"))
utils.mkifne(plotsdir("ripples","mpp_decode","withBehVideo=$usevideo"))
wells = task[(task.name.=="welllocs") .& (task.epoch .== epoch), :]


utils.pushover("Finished preprocess sequenece")

# -----------
# GET DECODER
# -----------
decode_file = raw.decodepath(animal, day, epoch, transition="empirical",
                             method="sortedspike", split=split,
                             type=split_type, speedup=20.0)
@assert isfile(decode_file) "File=$decode_file doesn't exist"
D = raw.load_decode(decode_file)
x = D["x_position"]
y = D["y_position"]
T = D["time"]
if usevideo
    stream = VideoIO.open(video); vid = VideoIO.openvideo(stream);
end
dat = permutedims(D[variable], [2,1,3])
#mint = minimum(beh[beh.epoch.==epoch,:].time)

# Ready theta waves

@showprogress for (cyc, cycle) in enumerate(eachrow(validcycles))
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
    ax = Axis(fig[4,1])
    α = 0.1
    for (i,unit) in enumerate(sp)
        try
            cmap = get(ColorSchemes.hawaii, ((unit.time.-start)./(stop-start)))
            cmap = parse.(Colorant, cmap)
            jitter= (α * rand(Float64,size(unit.time))) .- α/2
            scatter!(ax, unit.time, i*ones(size(unit.time))*4, color=cmap, strokewidth=1, markersize=3)
        catch
        end
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
    P(t) = Point2f.(behx(t), behy(t))
    area = ripple.area
    color         = decode.scatter.marker_color([ripple.area])
    sc_glow_color = decode.scatter.glow_color([ripple.area])
    sc_glow_width = ripple.amp

    # ----------------------------
    # SETUP FIGURE, AXIS, ARTISTIS
    # ----------------------------
    ax = Axis(fig[2:3,1], xlabel="centimeter", ylabel="centimeter", title="area=$(ripple.area), amp=$(round(ripple.amp,digits=2))")
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
        for e in ["pdf","svg"]
            savefile = plotsdir("ripples","mpp_decode", "withBehVideo=$usevideo",
                                 outputVideo,
                                 "rip=$rip.$area.amp=$(round(ripple.amp,digits=2))" *
                                 ".$e")
            save(savefile, fig, pt_per_unit = 1)
        end
    end

end

