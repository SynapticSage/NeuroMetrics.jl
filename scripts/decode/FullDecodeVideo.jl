using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
using Pushover, Revise, Interact, Blink, Mux, ProgressMeter
using Statistics
using VideoIO
using ColorSchemes
using GLMakie, Makie
using ColorSchemes, Colors
import ColorSchemeTools
set_theme!(theme_dark())
__revise_mode__ = :eval
#includet(srcdir("table.jl"))
includet(srcdir("raw.jl"))
includet(srcdir("table.jl"))
includet(srcdir("raster.jl"))
includet(srcdir("decode.jl"))
includet("/home/ryoung/Code/projects/goal-code/src/utils/SearchSortedNearest.jl/src/SearchSortedNearest.jl")
# -----------
# DATASETS
# -----------
animal, day, epoch, ca1tetrode = "RY16", 36, 7, 5
@time beh    = raw.load_behavior(animal, day)
@time spikes = raw.load_spikes(animal,   day)
@time cells  = raw.load_cells(animal,    day)
@time lfp    = raw.load_lfp(animal,      day)
@time task   = raw.load_task(animal,     day)
@time ripples= raw.load_ripples(animal,  day)
@time cycles = raw.load_cycles(animal,   day, tetrode)
wells = task[(task.name.=="welllocs") .& (task.epoch .== epoch), :]
video="/Volumes/Colliculus/RY16_experiment/actualVideos/RY16_69_0$(epoch)_CMt.1.mp4"

# -----------
# SETTINGS
# -----------
usevideo=false
dir = plotsdir("mpp_decode", "withBehVideo=$usevideo")
if !(isdir(dir))
    mkdir(dir)
end
transition_type, decoder_type = "empirical", "sortedspike"
variable                      = "causal_posterior"
outputVideo = "animation.$(decoder_type)_$(transition_type)_$(split_type)_$(split)_$(variable)_$(basename(video))"
(split, split_type) = collect(Iterators.product([0,1,2],
                                                ["testing","training"]))[1]
thresh = Dict("likelihood"=>0.1,
              "acausal_posterior"=>0.65,
              "causal_posterior"=> 0.65)

# -----------
# GET DECODER
# -----------
#D = raw.load_decode("/Volumes/FastData/RY16deepinsight36-07.$split_type.$transition_type.speedup=1.0.$split.nc")
#D = raw.load_decode("/Volumes/FastData/RY16deepinsight36-07.prediction.spikedecoder.test.nc")
decode_file = raw.decodepath(animal, day, epoch, transition="empirical",
                             method="sortedspike", split=split,
                             type=split_type, speedup=20.0)
D = raw.load_decode(decode_file)
beh, spikes, ripples, D = raw.normalize_time(beh, spikes, ripples, D; timefields=Dict(3=>["start", "stop", "time"]));
x = D["x_position"]
y = D["y_position"]
T = D["time"]
dat = D[variable]
D = nothing
stream = VideoIO.open(video)
vid = VideoIO.openvideo(stream)

# -------------------------------------------------
# Create a merged table of ripples and theta cycles
# -------------------------------------------------
events = vcat(transform(cycles,:end=>:stop), ripples, cols=:union, source=:source=>[:theta, :ripple]);
events = sort(events, :start)[!, Not(:end)]


# -----------------------------------------------
# Setup observable and connect decoder to behavior
# -----------------------------------------------
decode_inds = SearchSortedNearest.searchsortednearest.([beh.time], T)
t = Observable{Int}(1)

mint = minimum(beh[beh.epoch.==epoch,:].time)
dat = permutedims(dat, [2,1,3])
dat = decode.quantile_threshold(dat, thresh[variable])

#@time tetrode = raw.load_tetrode(animal,   day)
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
P             = @lift(Point2f.(behx($t), behy($t)))
color         = @lift(decode.scatter.color($t))
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
