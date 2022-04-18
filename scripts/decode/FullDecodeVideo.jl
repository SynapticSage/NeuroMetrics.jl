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
@time beh    = raw.load_behavior(animal, day)
@time spikes = raw.load_spikes(animal,   day)
@time cells  = raw.load_cells(animal,    day)
@time lfp    = @subset(raw.load_lfp(animal, day), :tetrode.==ca1_tetrode)[!,[:time, :raw, :phase, :amp]]
@time task   = raw.load_task(animal,     day)
@time ripples= raw.load_ripples(animal,  day)
@time cycles = raw.load_cycles(animal,   day, ca1_tetrode)
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
beh, lfp, spikes, ripples, D = raw.normalize_time(beh, lfp, spikes, ripples, D; 
                                                  timefields=Dict(4=>["start", "stop", "time"]));
x = D["x_position"]
y = D["y_position"]
T = D["time"]
dat = D[variable]
dat = permutedims(dat, [2,1,3])
dat = decode.quantile_threshold(dat, thresh[variable])
D = nothing
stream = VideoIO.open(video)
vid = VideoIO.openvideo(stream)
mint = minimum(beh[beh.epoch.==epoch,:].time)

# -------------------------------------------------
# Create a merged table of ripples and theta cycles
# -------------------------------------------------
lfp = begin
    lfp = raw.lfp.annotate_cycles(lfp, method="peak-to-peak")
    lfp.phase = raw.lfp.phase_to_radians(lfp.phase)
    beh, lfp = raw.register(beh, lfp; transfer=["velVec"], on="time")
    cycles = raw.lfp.get_cycle_table(lfp, :velVec=>mean)
    cycles = filter(:amp_mean => amp->(amp .> 50) .& (amp .< 600), cycles)
    cycles = filter(:δ => dur->(dur .> 0.025) .& (dur .< 0.4), cycles)
    cycles = filter(:velVec_mean => (𝒱  -> abs.(𝒱)  .> 2) , cycles)

    # Throw away bad cycles
    lfpcyc, cyc = lfp.cycle, cycles.cycle
    @time LoopVectorization.@avxt goodcycles = in.(lfpcyc, [cyc])
    lfp[(!).(goodcycles), :cycle] .= -1
    # Annotate rippples
    lfp.cycle, lfp.phase, lfp.time = Int32.(lfp.cycle), Float32.(lfp.phase), Float32.(lfp.time)
    ripples.type = ripples.area .* " ripple"
    ripples.rip_id = 1:size(ripples,1)
    lfp = raw.registerEvents(ripples, lfp, on="time", 
                             eventStart="start", eventStop="stop", 
                             transfer=["rip_id","type"])
    lfp.raw = Float32.(utils.norm_extrema(lfp.raw, extrema(spikes.unit)))
    lfp
end


t = Observable(Int(1))
Δt = median(diff(T))
Δ_bounds = [0.10, 0.10] # seconds
Δ = -Int(round(Δ_bounds[1]/Δt)) : Int(round(Δ_bounds[2]/Δt))

function select_range(t, data=spikes, Δ_bounds=Δ_bounds)
    time  = T[t]
    data = @subset(data, :time .> (time - Δ_bounds[1]),
                         :time .< (time + Δ_bounds[2]))
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

# Figure
Fig = Figure()
# Grids
gNeural = Fig[2,1] = GridLayout(1,1)
# Axes
@time behavior     = @lift select_time($t, beh)
@time title_str = @lift "corr=$($behavior.correct), goal=$($behavior.stopWell)"
axArena  = Axis(Fig[1,1], xlabel="x", ylabel="y", title=title_str)
axNeural = Axis(gNeural[1,1], xlabel="time")
@time spike_events = @lift select_range($t, spikes)
@time lfp_events   = @lift select_range($t, lfp)

# Plot elements

sc_sp_events = @lift([Point2f(Tuple(x)) for x in
                      eachrow($spike_events[!,[:time,:unit]])])
sc = scatter!(axNeural, sc_sp_events)

ln_lfp_events = @lift([Point2f(Tuple(x)) for x in
                       eachrow($lfp_events[!,[:time,:raw]])])
sc = lines!(axNeural, ln_lfp_events, color=:white)

hm_prob = @lift select_prob($t)
hm = heatmap!(axArena, x,y, hm_prob)
point_xy = @lift(Point2f($behavior.x, $behavior.y))
hm = scatter!(axArena, point_xy, color=:white)


# -------------
# CREATE VIDEO
# -------------
framerate = 60
timestamps = range(1, length(T), step=1)
#recording = plotsdir("mpp_decode", "withBehVideo=$usevideo", outputVideo)
record(Fig, "test.mp4", timestamps; framerate=framerate) do stamp
    t[] = stamp
end
