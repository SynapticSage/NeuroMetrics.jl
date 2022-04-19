using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
using Pushover, Revise, Interact, Blink, Mux, ProgressMeter
using Statistics
using VideoIO
using ColorSchemes
using GLMakie
using ColorSchemes, Colors
using DataFrames, DataFramesMeta
import ColorSchemeTools, LoopVectorization
using Printf
set_theme!(theme_dark())
__revise_mode__ = :eval
#includet(srcdir("table.jl"))
includet(srcdir("raw.jl"))
includet(srcdir("table.jl"))
includet(srcdir("raster.jl"))
includet(srcdir("decode.jl"))
includet(srcdir("utils.jl"))
includet("/home/ryoung/Code/projects/goal-code/src/utils/SearchSortedNearest.jl/src/SearchSortedNearest.jl")

# -----------------
# HELPER FUNCTIONS 
# -----------------
na = [CartesianIndex()]

# -----------
# DATASETS
# -----------
animal, day, epoch, ca1_tetrode = "RY16", 36, 7, 5
thresh = Dict("likelihood"=>0.1, "acausal_posterior"=>0.90, "causal_posterior"=> 0.90)
transition_type, decoder_type = "empirical", "sortedspike"
variable                      = "causal_posterior"
usevideo                      = false
remove_nonoverlap             = false # set to true if you expect multiple splits in this code
@time beh    = raw.load_behavior(animal, day)
@time spikes = raw.load_spikes(animal,   day)
@time cells  = raw.load_cells(animal,    day)
@time lfp    = @subset(raw.load_lfp(animal, day), :tetrode.==ca1_tetrode)[!,[:time, :raw, :phase, :amp]]
@time task   = raw.load_task(animal,     day)
@time ripples= raw.load_ripples(animal,  day)
wells = task[(task.name.=="welllocs") .& (task.epoch .== epoch), :]
video="/Volumes/Colliculus/RY16_experiment/actualVideos/RY16_69_0$(epoch)_CMt.1.mp4"

# -----------
# PREPROCESS
# -----------
dir = plotsdir("mpp_decode", "withBehVideo=$usevideo")
if !(isdir(dir))
    mkdir(dir)
end
outputVideo = "animation.$(decoder_type)_$(transition_type)_$(split_type)_$(split)_$(variable)_$(basename(video))"
(split, split_type) = collect(Iterators.product([0,1,2], 
                                                ["test","train"]))[1]

# -----------
# GET DECODER
# -----------
decode_file = raw.decodepath(animal, day, epoch, transition="empirical",
                             method="sortedspike", split=split,
                             type=split_type, speedup=20.0)
D = raw.load_decode(decode_file)
stream = VideoIO.open(video)
vid = VideoIO.openvideo(stream)

# Normalize times
mint = minimum(beh[beh.epoch.==epoch,:].time)
beh, spikes, lfp, ripples, D = raw.normalize_time(beh, spikes, lfp, ripples, D; 
                                                  timefields=Dict(4=>["start", "stop", "time"]));
x, y, T, dat = D["x_position"], D["y_position"], D["time"], D[variable]
dat = permutedims(dat, [2,1,3])
dat = decode.quantile_threshold(dat, thresh[variable])
D = nothing
nmint = minimum(beh[beh.epoch.==epoch,:].time)

# -------------------------------------------------
# Create a merged table of ripples and theta cycles
# -------------------------------------------------
# (1) Annotate and filter Θ cycles 
lfp = raw.lfp.annotate_cycles(lfp, method="peak-to-peak")
lfp.phase = raw.lfp.phase_to_radians(lfp.phase)
beh, lfp = raw.register(beh, lfp; transfer=["velVec"], on="time")
lfp.raw = Float32.(utils.norm_extrema(lfp.raw, extrema(spikes.unit)))

# (2) Throw away bad Θ cycles
cycles = raw.lfp.get_cycle_table(lfp, :velVec=>mean)
cycles = filter(:amp_mean => amp->(amp .> 50) .& (amp .< 600), cycles)
cycles = filter(:δ => dur->(dur .> 0.025) .& (dur .< 0.4), cycles)
cycles = filter(:velVec_mean => (𝒱  -> abs.(𝒱)  .> 2) , cycles)

# (3) Remove cycle labels in lfp of bad Θ cycles
#lfpcyc, cyc = lfp.cycle, cycles.cycle
lfp.chunks = Int16.(round.(((1:size(lfp,1))./1000000)))
lfp = groupby(lfp,:chunks)
@time Threads.@threads for group in lfp
    LoopVectorization.@avxt goodcycles = in.(group.cycle, [cycles.cycle])
    group[(!).(goodcycles), :cycle] .= -1
end
lfp = combine(lfp, identity)

# (4) Annotate ripples into lfp
lfp.cycle, lfp.phase, lfp.time = Int32.(lfp.cycle), Float32.(lfp.phase),
                                 Float32.(lfp.time)
ripples.type = ripples.area .* " ripple"
ripples.rip_id = 1:size(ripples,1)
lfp = raw.registerEvents(ripples, lfp, 
                         on="time", 
                         eventStart="start", 
                         eventStop="stop", 
                         transfer=["rip_id","type"])

if remove_nonoverlap
    spikes, beh, lfp, ripples, T_inds = raw.keep_overlapping_times(spikes, beh, lfp, ripples, T;
                                                                  timefields=Dict(4=>["start", "stop", "time"],
                                                                  returninds=[5])
    dat, T = dat[:,:,T_inds], T[T_inds]
end


t = Observable(Int(1))
Δt = median(diff(T))
Δ_bounds = [0.20, 0.20] # seconds
Δ = -Int(round(Δ_bounds[1]/Δt)) : Int(round(Δ_bounds[2]/Δt))

function select_range(t, data=spikes, Δ_bounds=Δ_bounds)
    time  = T[t]
    data = @subset(data, (:time .> (time - Δ_bounds[1])) .&&
                           (:time .< (time + Δ_bounds[2])))
    data.time .-= T[t]
    data
end
function select_time(t, data=spikes, Δ_bounds=Δ_bounds)
    time  = T[t]
    I = utils.searchsortednearest(beh.time, time)
    data = data[I, :]
    data.time = 0
    data
end
function select_prob(t, prob=dat)
    time  = min(max(t, 1), length(T))
    D = prob[:,:, time]
end
@time behavior     = @lift select_range($t, beh, [0.20, 0.01])
@time spike_events = @lift select_range($t, spikes)
@time lfp_events   = @lift select_range($t, lfp)

# Figure
Fig = Figure()
# Grids
gNeural = Fig[2,1] = GridLayout(1,1)
# Axes
@time title_str = @lift "corr=$($behavior.correct[1]), goal=$($behavior.stopWell[1]), vel=$(@sprintf("%2.1f",abs($behavior.velVec[1])))"
axArena  = Axis(Fig[1,1], xlabel="x", ylabel="y", title=title_str)
axNeural = Axis(gNeural[1,1], xlabel="time")
# Plot elements
sc_sp_events = @lift([Point2f(Tuple(x)) for x in
                      eachrow($spike_events[!,[:time,:unit]])])
sc = scatter!(axNeural, sc_sp_events, color=:white, markersize=3)
ln_lfp_events = @lift([Point2f(Tuple(x)) for x in
                       eachrow($lfp_events[!,[:time,:raw]])])
ln_lfp_phase = @lift([Point2f(Tuple(x)) for x in
                      eachrow(select($lfp_events, :time, :phase′=>x->x.*30))])
sc = lines!(axNeural, ln_lfp_events, color=:white)
sc = lines!(axNeural, ln_lfp_phase,  color=:white, linestyle=:dash)
hm_prob = @lift select_prob($t)
hm = heatmap!(axArena, x,y, hm_prob, colormap=(:bamako,0.8))
point_xy = @lift(Point2f($behavior.x[1], $behavior.y[1]))
line_xy = @lift([Point2f(b.x, b.y) for b in eachrow($behavior)])
sc_xy = scatter!(axArena, point_xy, color=:white)
ln_xy = lines!(axArena, line_xy, color=:white, linestyle=:dash)


# -------------
# CREATE VIDEO
# -------------
framerate = 60
timestamps = range(100000, length(T), step=4)
#recording = plotsdir("mpp_decode", "withBehVideo=$usevideo", outputVideo)
record(Fig, "test.mp4", timestamps; framerate=framerate) do stamp
    t[] = stamp
end
