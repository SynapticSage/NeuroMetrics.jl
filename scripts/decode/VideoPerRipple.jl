using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
using Pushover, Revise
using Interact, Blink, Mux, ProgressMeter
using ProgressMeter
using Statistics
using VideoIO
using ColorSchemes
using GLMakie, Makie
set_theme!(theme_dark())
__revise_mode__ = :eval
#includet(srcdir("table.jl"))
includet(srcdir("raw.jl"))
includet(srcdir("table.jl"))
includet(srcdir("raster.jl"))
includet(srcdir("decode.jl"))
includet("/home/ryoung/Code/projects/goal-code/src/utils/SearchSortedNearest.jl/src/SearchSortedNearest.jl")
beh     = raw.load_behavior("RY16", 36)
spikes  = raw.load_spikes("RY16",   36)
ripples = raw.load_ripples("RY16",  36)
cmtopx(x) = 0.1487 * x
pxtocm(x) = x / 0.1487 

# -----------
# SETTINGS
# -----------
usevideo=false
mkdir(plotsdir("ripples","mpp_decode"))
mkdir(plotsdir("ripples","mpp_decode","withBehVideo=$usevideo"))
transition_type ="empirical"
decoder_type ="sortedspike"
variable = "causal_posterior"
epoch = 7
video="/Volumes/Colliculus/RY16_experiment/actualVideos/RY16_69_0$(epoch)_CMt.1.mp4"
(split, split_type) = collect(Iterators.product([0,1,2], ["testing","training"]))[1]
thresh = Dict("likelihood"=>0.1,
              "acausal_posterior"=>0.65,
              "causal_posterior"=>0.65)
outputVideo = "animation.$(decoder_type)_$(transition_type)_$(split_type)_$(split)_$(variable)_$(basename(video))"
mkdir(plotsdir("ripples","mpp_decode", "withBehVideo=$usevideo", outputVideo))

# -----------
# GET DECODER
# -----------
#decoder = "/Volumes/FastData/RY16deepinsight36-07.$split_type.$transition_type.speedup=1.0.$split.nc"
#decode_file = "/Volumes/FastData/RY16deepinsight36-07.prediction.spikedecoder.test.nc"
#decode_file = "/Volumes/FastData/decode_unbinned/RY16deepinsight36-0$epoch.$split_type.$decoder_type.$transition_type.speedup=1.0.$split.nc"
#decode_file = "/Volumes/FastData/decode_unbinned/RY16deepinsight36-0$epoch.$split_type.$decoder_type.$transition_type.speedup=1.0.$split.nc"
decode_file = "/mnt/FastData/RY16deepinsight36-07.testing.sortedspike.empirical.notbinned.downsamp=1.speedup=20.0.iter=0.nc"
D = raw.load_decode(decode_file)
x = D["x_position"]
y = D["y_position"]
T = D["time"]
if usevideo
    stream = VideoIO.open(video)
    vid = VideoIO.openvideo(stream)
end

mint = minimum(beh[beh.epoch.==epoch,:].time)
dat = permutedims(D[variable], [2,1,3])
#dat = decode.quantile_threshold(dat, thresh=thresh[variable])

validripples = ripples[ripples.epoch.==epoch,:]
#(rip,ripple) = collect(enumerate(eachrow(validripples)))[4]
for (rip,ripple) in enumerate(eachrow(validripples))

    start, stop = ripple.start, ripple.stop
    B = beh[(beh.time.>=start) .& (beh.time.<=stop),:]
    dat_sub = dat[:,:,(D["time"].>=start) .& (D["time"].<=stop)]
    if size(dat_sub,3) == 0
        continue
    end

    # -----------------------------------------------
    # Setup observable and connect decoder to behavior
    # -----------------------------------------------
    t = Observable{Int}(1)
    decode_inds = SearchSortedNearest.searchsortednearest.([B.time], D["time"])

    # -----------------
    # READY HEATMAP VARS
    # -----------------
    H(t) = dat_sub[:,:,t]
    HH = @lift(H($t))

    # -----------------
    # READ SCATTER VARS
    # -----------------
    behx(t) = B[decode_inds[t],"x"]
    behy(t) = B[decode_inds[t],"y"]
    behpoint(t) = (behx(t), behy(t))
    P = @lift(Point2f.(behx($t), behy($t)))
    area = ripple.area
    color         = decode.scatter.marker_color([ripple.area])
    sc_glow_color = decode.scatter.glow_color([ripple.area])
    sc_glow_width = ripple.amp

    # ----------------------------
    # SETUP FIGURE, AXIS, ARTISTIS
    # ----------------------------
    if usevideo
        FRAME = @lift(decode.video.frameattime(vid, T[$t]-mint-(0*0.033)))
        fig=image(x,y,FRAME, depth_shift=0)
        hm_=heatmap!(x,y,HH, transparency=true, depth_shift=0.5,
                     colormap=(:bamako, 0.25))
    else
        fig = Figure()
        ax = Axis(fig[1,1], xlabel="cm", ylabel="cm", title="area=$(ripple.area), amp=$(round(ripple.amp,digits=2))")
        FRAME = nothing
        hm_=heatmap!(ax, x,y,HH, transparency=true, depth_shift=0.5,
                     colormap=(:bamako, 0.5))
    end
    sc_ = scatter!(P, color=color, depth_shift=0.8, overdraw=true,
                   markersize=20, glowwidth=sc_glow_width, glowcolor=(sc_glow_color, 0.8),
                   label="position")
    if usevideo
        Legend(fig.figure[1, 2], [sc_], ["Actual\nposition"])
    else
        Legend(fig[1, 2], [sc_], ["Actual\nposition"])
    end

    # -------------
    # CREATE VIDEO
    # -------------
    framerate = 60
    timestamps = range(1, size(dat_sub,3), step=1)
    recording = plotsdir("ripples","mpp_decode", "withBehVideo=$usevideo", outputVideo, "rip=$rip.$area.amp=$(round(ripple.amp,digits=2))" * ".mp4")
    record(fig, recording, timestamps; framerate=framerate) do stamp
        t[] = stamp 
    end

end
