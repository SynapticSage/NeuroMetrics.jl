using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
using Pushover, Revise, Interact, Blink, Mux, ProgressMeter
using Statistics, NaNStatistics
using VideoIO
using ColorSchemes
using GLMakie
using ColorSchemes, Colors
using DataFrames, DataFramesMeta
import ColorSchemeTools, LoopVectorization
using Printf
using StatsPlots: @df
set_theme!(theme_dark())
__revise_mode__ = :eval
#includet(srcdir("table.jl"))
includet(srcdir("raw.jl"))
includet(srcdir("table.jl"))
includet(srcdir("raster.jl"))
includet(srcdir("decode.jl"))
includet(srcdir("utils.jl"))
import raw, table, raster, decode, utils
import StatsBase
ENV["JULIA_DEBUG"] = nothing

# -----------------
# HELPER FUNCTIONS 
# -----------------
na = [CartesianIndex()]

# -----------
# DATASETS
# -----------
animal, day, epoch, ca1_tetrode = "RY16", 36, 7, 5
thresh = Dict("likelihood"=>0.1, "acausal_posterior"=>0.95, "causal_posterior"=> 0.95)
transition_type, decoder_type = "empirical", "sortedspike"
variable                      = "causal_posterior"
usevideo                      = false
remove_nonoverlap             = true # set to true if you expect multiple splits in this code
dosweep                       = false
@time beh    = raw.load_behavior(animal, day)
@time spikes = raw.load_spikes(animal,   day)
@time cells  = raw.load_cells(animal,    day)
@time lfp    = @subset(raw.load_lfp(animal, day), :tetrode.==ca1_tetrode)[!,[:time, :raw, :phase, :amp]]
@time task   = raw.load_task(animal,     day)
@time ripples= raw.load_ripples(animal,  day)
wells = task[(task.name.=="welllocs") .& (task.epoch .== epoch), :]
homeWell, arenaWells = begin
    hw = argmax(StatsBase.fit(StatsBase.Histogram, filter(b->b!=-1,beh.stopWell), 1:6).weights)
    wells[hw,:], wells[setdiff(1:5,hw),:]
end
boundary = task[(task.name.=="boundary") .& (task.epoch .== epoch), :]
append!(boundary, DataFrame(boundary[1,:]))
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
beh, ripples = raw.register(beh, ripples; transfer=["velVec"], on="time")
ripples = ripples[abs.(ripples.velVec) .< 2, :]
lfp.raw = Float32.(utils.norm_extrema(lfp.raw, extrema(spikes.unit)))

# (2) Throw away bad Θ cycles
cycles = raw.lfp.get_cycle_table(lfp, :velVec => (x->median(abs.(x))) => :velVec_median )
cycles = filter(:amp_mean => amp->(amp .> 50) .& (amp .< 600), cycles)
cycles = filter(:δ => dur->(dur .> 0.025) .& (dur .< 0.4), cycles)
cycles = filter(:velVec_median => (𝒱  -> abs.(𝒱)  .> 2) , cycles)
# TODO Remove any cycles in side a ripple

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
                                                                  returninds=[5])
    dat, T = dat[:,:,T_inds], T[T_inds]
end
[extrema(x.time) for x in (lfp, spikes, ripples, beh)]


# Add ripple phase
lfp = sort(combine(lfp, identity), :time)
lfp.rip_phase = Float32.(combine(groupby(lfp, :rip_id, sort=false), x->1/nrow(x)*ones(nrow(x))).x1)
# Theta : Create probability chunks by phase
dat = Float32.(dat)
theta, ripple, non = copy(dat), copy(dat), copy(dat)
@time Threads.@threads for (t,time) in collect(enumerate(T))
    I = utils.searchsortednearest(lfp.time, time)
    θ = theta[:, :, t] 
    ρ = ripple[:, :, t] 
    not_a_theta_cycle = lfp.cycle[I] == -1
    is_a_ripple = !(ismissing(lfp.rip_id))
    if not_a_theta_cycle # NOT THETA
        θ .= NaN
        if is_a_ripple # IS RIPPLE?
            ρ[(!).(isnan.(θ))] .= lfp.rip_phase[I]
            non[:,:,t] .= NaN
        end
    else # THETA CYCLE
        θ[(!).(isnan.(θ))] .= lfp.phase[I]
        ρ .= NaN
        non[:,:,t] .= NaN
    end
    theta[:,:,t] = θ
    ripple[:,:,t] = ρ
end

utils.pushover("Ready to plot")

if dosweep

    # Create cumulative theta sweeps
    sweep = (a,b)->isnan(b) ? a : nanmean(cat(a, b, dims=3), dims=3)
    lfp = groupby(lfp,:cycle)
    @time @Threads.threads for group in lfp
        cycStart, cycStop = utils.searchsortednext(T, group.time[1]),
                            utils.searchsortednext(T, group.time[end])
        local cycle = cycStart:cycStop
        if cycle == 1:1 || group.cycle[1] == -1 || ismissing(group.cycle[1])
            continue
        end
        theta[:, :, cycle] = accumulate((a,b)->sweep.(a,b), theta[:,:,cycle],  dims=3)
    end
    lfp = combine(lfp, identity)
    lfp = sort!(lfp,:time)

    # Cumulative ripple sweeps
    lfp = groupby(lfp,:rip_id)
    @time @Threads.threads for group in lfp
        cycStart, cycStop = utils.searchsortednext(T, group.time[1]), utils.searchsortednext(T, group.time[end])
        local cycle = cycStart:cycStop
        if cycle == 1:1 || group.cycle[1] == -1 || ismissing(group.cycle[1])
            continue
        end
        ripple[:, :, cycle] = accumulate((a,b)->sweep.(a,b), ripple[:,:,cycle],  dims=3)
    end
    lfp = combine(lfp,identity)

end



Δ_bounds = [0.20, 0.20] # seconds
tr = Dict("beh"=>utils.searchsortednearest(beh.time, T[1]),
          "lfp"=>utils.searchsortednearest(beh.time, T[1]))
Δt = Dict("beh"=>median(diff(beh.time)),
          "lfp"=>median(diff(lfp.time)),
          "prob"=>median(diff(T)))
Δi = Dict("beh"=>(Δt["beh"]/Δt["prob"])^-1,
          "lfp"=>(Δt["lfp"]/Δt["prob"])^-1)

function select_range(t, data=spikes, Δ_bounds=Δ_bounds)
    time  = T[t]
    data = @subset(data,   (:time .> (time - Δ_bounds[1])) .&&
                           (:time .< (time + Δ_bounds[2])))
    data.time = data.time .- T[t]
    data
end
function select_est_range(t, tr, Δt, Δi, data=beh, Δ_bounds=Δ_bounds)
    tt = tr + (t-1)*Δi
    Δ = -Int(round(Δ_bounds[1]/Δt)) : Int(round(Δ_bounds[2]/Δt))
    center_time = data.time[t]
    data = data[tt+Δ[1]:tt+Δ[2],:]
    data.time .-= center_time
    data
end
function select_est_range(t, tr, Δt, data=beh, Δ_bounds=Δ_bounds)
    I = utils.searchsortednearest(data.time, T[t])
    if I != 1 && !(isnan(I))
        Δ = -Int(round(Δ_bounds[1]/Δt)) : Int(round(Δ_bounds[2]/Δt))
        center_time = data.time[I]
        @debug "I=$I, Δ=$Δ"
        data = data[max.(I.+Δ,1),:]
        data.time .-= center_time
    else
        data = DataFrame(lfp[1,:])
        data.time .= NaN
    end
    return data
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

## FIGURE WISHLIST
# - Goals
#   - Glow goal
#   - Vector to goal?
# - Sequences
#   - Start to end
#   - Scatter of maxima
# - Boundaries
# - Home well

## FIGURE SETTINGS
visualize = :slider

# Figure
Fig = Figure()
if visualize == :video
    t = Observable(Int(1000))
elseif visualize == :slider
    sl_t = Slider(Fig[1:2, 2], range = 1:length(T), horizontal = false, 
                  startvalue = 1000)
    t = lift(sl_t.value) do tt
        tt
    end
end

# Data
@time behavior     = @lift select_range($t, beh, [-0.01, 0.60])
#@time behavior     = @lift select_est_range($t, tr["beh"], Δt["beh"], beh, [0.60, 0.01])
@time spike_events = @lift select_range($t, spikes)
#@time lfp_events   = @lift select_range($t, lfp)
@time lfp_events   = @lift select_est_range($t, tr["lfp"], Δt["lfp"], lfp)

# Grids
gNeural = Fig[2,1] = GridLayout(1,1)
controls = Fig[1:2, 3]
# Axes
@time title_str = @lift "corr=$($behavior.correct[1]), goal=$($behavior.stopWell[1]), vel=$(@sprintf("%2.1f",abs($behavior.velVec[1])))"
axArena  = Axis(Fig[1,1], xlabel="x", ylabel="y", title=title_str)
axNeural = Axis(gNeural[1,1], xlabel="time")

# Neural data AXIS
sc_sp_events = @lift([Point2f(Tuple(x)) for x in
                      eachrow($spike_events[!,[:time,:unit]])])
sc = scatter!(axNeural, sc_sp_events, color=:white, markersize=3)
ln_lfp_events = @lift([Point2f(Tuple(x)) for x in eachrow($lfp_events[!,[:time,:raw]])]) sc = lines!(axNeural, ln_lfp_events, color=:white)
ln_lfp_phase = @lift([Point2f(Tuple(x)) for x in
                      eachrow(select($lfp_events, :time, :phase=>x->x.*30))])
sc = lines!(axNeural, ln_lfp_phase,  color=:gray, linestyle=:dash)
vlines!(axNeural, [0], color=:red, linestyle=:dash)

# Arena data AXIS
#hm_prob = @lift select_prob($t)
#hm = heatmap!(axArena, x,y, hm_prob, colormap=(:bamako,0.8))
theta_prob = @lift select_prob($t, theta)
hm_theta = heatmap!(axArena, x,y, theta_prob, colormap=(:buda,0.8), colorrange=(-pi, pi), interpolate=false)
ripple_prob = @lift select_prob($t, ripple)
hm_ripple = heatmap!(axArena, x,y, ripple_prob, colormap=(:diverging_tritanopic_cwr_75_98_c20_n256, 0.8), interpolate=false)
#non_prob = @lift select_prob($t, non)
#hm_non = heatmap!(axArena, x,y, non_prob, colormap=(:bamako,0.8))
now = 1
point_xy = @lift(Point2f($behavior.x[now], $behavior.y[now]))
line_xy = @lift([Point2f(b.x, b.y) for b in eachrow($behavior)])
sc_xy = scatter!(axArena, point_xy, color=:white)
sc_arena = scatter!(axArena, arenaWells.x, arenaWells.y, marker=:star5, markersize=40, color=:gray)
sc_home = scatter!(axArena, [homeWell.x], [homeWell.y],   marker='𝐇', markersize=25,   color=:gray)
future = @lift $behavior.stopWell[now]
future_well_xy = @lift $future == -1 ? Point2f(NaN, NaN) : Point2f(wells.x[$future], wells.y[$future])
sc_future_well = scatter!(axArena, future_well_xy, marker='𝐅', markersize=25, color=:white, glowwidth=5)
#past = @lift $behavior.stopWell[now]
sc_future = scatter!(axArena, [homeWell.x], [homeWell.y], marker='𝐇', markersize=25, color=:white)
ln_xy = lines!(axArena, line_xy, color=:white, linestyle=:dash)
lines!(axArena, boundary.x, boundary.y, color=:grey)

function play()
    framerate = 60
    timestamps = range(t[], length(T), step=1)
    #recording = plotsdir("mpp_decode", "withBehVideo=$usevideo", outputVideo)
    record(Fig, "test.mp4", timestamps; framerate=framerate) do stamp
        t[] = stamp
    end
end


# -------------
# VISUALIZE
# -------------
if visualize == :video
    play()
elseif visualize == :slider
end

