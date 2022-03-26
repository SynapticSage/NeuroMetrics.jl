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
beh = raw.load_behavior("RY16",36)
ripples = raw.load_ripples("RY16",36)
cmtopx(x) = 0.1487 * x
pxtocm(x) = x/0.1487 

# -----------
# SETTINGS
# -----------
usevideo=false
mkdir(plotsdir("mpp_decode", "withBehVideo=$usevideo"))
transition_type ="empirical"
decoder_type ="sortedspike"
variable = "causal_posterior"
epoch = 7
video="/Volumes/Colliculus/RY16_experiment/actualVideos/RY16_69_0$epoch_CMt.1.mp4"
outputVideo = "animation.$(decoder_type)_$(transition_type)_$(split_type)_$(split)_$(variable)_$(basename(video))"
(split, split_type) = collect(Iterators.product([0,1,2], ["testing","training"]))[1]
thresh = Dict("likelihood"=>0.1,
              "acausal_posterior"=>0.65,
              "causal_posterior"=>0.65)

# -----------
# GET DECODER
# -----------
#D = raw.load_decode("/Volumes/FastData/RY16deepinsight36-07.$split_type.$transition_type.speedup=1.0.$split.nc")
#D = raw.load_decode("/Volumes/FastData/RY16deepinsight36-07.prediction.spikedecoder.test.nc")
D = raw.load_decode("/Volumes/FastData/RY16deepinsight36-0$epoch.$split_type.$decoder_type.$transition_type.speedup=1.0.$split.nc")
x = D["x_position"]
y = D["y_position"]
T = D["time"]
stream = VideoIO.open(video)
vid = VideoIO.openvideo(stream)

# -----------------------------------------------
# Setup observable and connect decoder to behavior
# -----------------------------------------------
decode_inds = SearchSortedNearest.searchsortednearest.([beh.time], D["time"])
t = Observable{Int}(1)

mint = minimum(beh[beh.epoch.==epoch,:].time)
dat = permutedims(D[variable], [2,1,3])
dat = decodee.quantile_threshold(dat)

# -----------------
# READY HEATMAP VARS
# -----------------
H(t) = dat[:,:,t]
HH = @lift(H($t))

# -----------------
# READ SCATTER VARS
# -----------------
behx(t) = beh[decode_inds[t],"x"]
behy(t) = beh[decode_inds[t],"y"]
behpoint(t) = (behx(t), behy(t))
P = @lift(Point2f.(behx($t), behy($t)))
color = @lift(decoode.scatter.color($t))
sc_glow_color = @lift(decode.scatter.glow_color($t))
sc_glow_width = @lift(decode.scatter.glow_width($t))

# ----------------------------
# SETUP FIGURE, AXIS, ARTISTIS
# ----------------------------
if usevideo
    FRAME = @lift(raw.video.frameattime(vid, T[$t]-mint-(0*0.033)))
    fig=image(x,y,FRAME, depth_shift=0)
    hm_=heatmap!(x,y,HH, transparency=true, depth_shift=0.5,
                 colormap=(:bamako, 0.25))
else
    fig = Figure()
    ax = Axis(fig[1,1], xlabel="cm", ylabel="cm", title="")
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
timestamps = range(1, length(T), step=1)
recording = plotsdir("mpp_decode", "withBehVideo=$usevideo", outputVideo)
record(fig, recording, timestamps; framerate=framerate) do stamp
    t[] = stamp 
end
