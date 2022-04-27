# -------
# IMPORTS
# -------
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
using Infiltrate
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
thresh = Dict("likelihood"=>0.1, "acausal_posterior"=>0.97, "causal_posterior"=> 0.97)
typical_number_of_squares_active = quantile(1:(45*28), 1) - quantile(1:(45*28), 0.97)
transition_type, decoder_type = "empirical", "sortedspike"
variable                      = "causal_posterior"
usevideo                      = false
remove_nonoverlap             = true # set to true if you expect multiple splits in this code
dosweep                       = false
doPrevPast                    = false
doRipplePhase                 = false
splitBehVar                   = ["egoVec_1", "egoVec_2", "egoVec_3", "egoVec_4", "egoVec_5"] # Variables to monitor with splits
vectorToWells                 = true
histVectorToWells             = true
@time beh    = raw.load_behavior(animal, day)
@time spikes = raw.load_spikes(animal,   day)
@time cells  = raw.load_cells(animal,    day)
@time lfp    = @subset(raw.load_lfp(animal, day), :tetrode.==ca1_tetrode)[!,[:time, :raw, :phase, :amp, :broadraw]]
@time task   = raw.load_task(animal,     day)
@time ripples= raw.load_ripples(animal,  day)
wells = task[(task.name.=="welllocs") .& (task.epoch .== epoch), :]
homeWell, arenaWells = begin
    hw = argmax(StatsBase.fit(StatsBase.Histogram, filter(b->b!=-1,beh.stopWell), 1:6).weights)
    @info "Homewell = $hw"
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
(split_num, split_type) = collect(Iterators.product([0,1,2], 
                                                ["test","train"]))[1]
outputVideo = "animation.$(decoder_type)_$(transition_type)_$(split_type)_$(split_num)_$(variable)_$(basename(video))"

# -----------
# GET DECODER
# -----------
decode_file = raw.decodepath(animal, day, epoch, transition="empirical",
                             method="sortedspike", split=split_num,
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
# PREPROCESS LFP 
# TODO Turn this into a function
# -------------------------------------------------
# (1) Annotate and filter Θ cycles 
lfp = raw.lfp.annotate_cycles(lfp, method="peak-to-peak")
lfp.phase = raw.lfp.phase_to_radians(lfp.phase)
beh, lfp = raw.register(beh, lfp; transfer=["velVec"], on="time")
beh, ripples = raw.register(beh, ripples; transfer=["velVec"], on="time")
ripples = ripples[abs.(ripples.velVec) .< 2, :]
lfp.raw = Float32.(utils.norm_extrema(lfp.raw, extrema(spikes.unit)))
lfp.broadraw = Float32.(utils.norm_extrema(lfp.broadraw, extrema(spikes.unit)))

# (2) Throw away bad Θ cycles
cycles = raw.lfp.get_cycle_table(lfp, :velVec => (x->median(abs.(x))) => :velVec_median)
transform!(cycles, [:start,:end] => ((x,y)->mean([x,y])) => :time)
cycles = filter(:amp_mean => amp->(amp .> 50) .& (amp .< 600), cycles) cycles = filter(:δ => dur->(dur .> 0.025) .& (dur .< 0.4), cycles)
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
lfp.cycle, lfp.phase = Int32.(lfp.cycle), Float32.(lfp.phase)
                                 
ripples.type = ripples.area .* " ripple"
ripples.rip_id = 1:size(ripples,1)
lfp = begin
    if "area" in names(lfp)
        lfp[!,Not("area")]
    else
        lfp
    end
end
lfp = raw.registerEventsToContinuous(ripples, lfp, 
                         on="time", 
                         eventStart="start", 
                         eventStop="stop", 
                         ifNonMissingAppend=true,
                         targetEltype=Dict("area"=>String),
                         transfer=["rip_id","type","area"])
lfp.phase_plot = utils.norm_extrema(lfp.phase, extrema(spikes.unit))
lfp.phase = utils.norm_extrema(lfp.phase, (-pi,pi))

# (5) Annotate cycles with, decode vector, next/previous goal data
match(time, col) = beh[utils.searchsortednearest(beh.time, time),col]
function matchdxy(time::Real) 
    @infiltrate
    I =  utils.searchsortednearest(T, time)
    D = replace(dat[:,:,I], NaN=>0)
    xi = argmax(maximum(D, dims=2), dims=1)
    yi = argmax(utils.squeeze(maximum(D, dims=1)), dims=1)
    [x[xi][1], y[yi][1]]
end

cycles = transform(cycles, :start=>(t->match.(t, "x"))=> :start_x,
                   :start=>(x->matchdxy.(x))=>[:start_x_dec, :start_y_dec],
                   :end=>(x->matchdxy.(x))=>[:end_x_dec, :end_y_dec])

if remove_nonoverlap
    spikes, beh, lfp, ripples, T_inds = raw.keep_overlapping_times(spikes, beh, lfp, ripples, T;
                                                                  returninds=[5])
    dat, T = dat[:,:,T_inds], T[T_inds]
end
[extrema(x.time) for x in (lfp, spikes, ripples, beh)]


# ------------------------------
# Separate DECODES by LFP event!
# TODO Turn this into a function
# ------------------------------
lfp = sort(combine(lfp, identity), :time)
lfp.rip_phase = Float32.(combine(groupby(lfp, :rip_id, sort=false),
                                 x->1/nrow(x)*ones(nrow(x))).x1)
# Theta : Create probability chunks by phase
dat = Float32.(dat)
theta, ripple = copy(dat), repeat(copy(dat), outer=(1,1,1,4))
@time Threads.@threads for (t,time) in collect(enumerate(T))
    I = utils.searchsortednearest(lfp.time, time)
    θ, ρ  = view(theta, :, :, t), view(ripple, :, :, t, :)
    not_a_theta_cycle = lfp.cycle[I] == -1
    is_a_ripple = !(ismissing(lfp.rip_id[I]))
    if not_a_theta_cycle # NOT THETA
        θ .= NaN
        if is_a_ripple # IS RIPPLE?
            if lfp.area[I] == "CA1"
                ρ[:, :,  [2, 3, 4]] .= NaN
            elseif lfp.area[I] == "CA1PFC"
                ρ[:, :,  [1, 3, 4]] .= NaN
            elseif lfp.area[I] == "PFCCA1"
                ρ[:, :,  [1, 2, 4]] .= NaN
            elseif lfp.area[I] == "PFC"
                ρ[:, :,  [1, 2, 3]] .= NaN
            else
                @error "Not a valid type"
            end

            if doRipplePhase
                ρ[(!).(isnan.(ρ))]   .= lfp.rip_phase[I]
            end
        else
            ρ .= NaN
        end
    else # THETA CYCLE
        θ[(!).(isnan.(θ))] .= lfp.phase[I]
        ρ .= NaN
    end
end
frac_print(frac) = @sprintf("Fraction times = %2.4f",frac)
frac = mean(all(isnan.(theta) .&& isnan.(ripple), dims=(1,2)))
frac_print(frac)
frac_rip = vec(mean(all(isnan.(ripple),dims=(1,2)), dims=(1,2,3)))
frac_print.(frac_rip)

if dosweep

    # Create cumulative theta sweeps
    sweep = (a,b)->isnan(b) ? a : nanmean(cat(a, b, dims=3), dims=3)
    lfp = groupby(lfp,:cycle)
    @time Threads.@threads for group in lfp
        cycStart, cycStop = utils.searchsortednext(T, group.time[1]),
                            utils.searchsortednext(T, group.time[end])
        cycle = cycStart:cycStop
        if cycle == 1:1 || group.cycle[1] == -1 || ismissing(group.cycle[1])
            continue
        end
        θ = view(theta, :, :, cycle)
        θ = accumulate(sweep, theta[:,:,cycle],  dims=3)
    end
    lfp = combine(lfp, identity)
    lfp = sort!(lfp,:time)

    # Cumulative ripple sweeps
    if doRipplePhase
        lfp = groupby(lfp, :rip_id)
        @time @Threads.threads for group in lfp
            cycStart, cycStop = utils.searchsortednext(T, group.time[1]),
                                utils.searchsortednext(T, group.time[end])
            local cycle = cycStart:cycStop
            if cycle == 1:1 || group.cycle[1] == -1 ||
               ismissing(group.cycle[1])
                continue
            end
            ripple[:, :, cycle] = accumulate((a,b)->sweep.(a,b),
                                             ripple[:,:,cycle], dims=3)
        end
        lfp = combine(lfp,identity)
    end

end

# --------------------
# Preprocess BEHAHVIOR
# TODO Turn this into a function
# --------------------
beh = groupby(beh, [:epoch, :traj])
for (g,group) in enumerate(beh)
    if g!=length(beh)
        group.futureStopWell .= beh[g+1].stopWell[1]
    end
    if g!=1
        group.pastStopWell .= beh[g-1].stopWell[1]
    end
end
beh = sort(combine(beh, identity), :time)
replace!(beh.pastStopWell, missing=>-1)
replace!(beh.futureStopWell, missing=>-1)
if doPrevPast
    beh.hw = argmax(StatsBase.fit(StatsBase.Histogram,
                                  filter(b->b!=-1,beh.stopWell), 1:6).weights)
    beh = groupby(beh, [:epoch, :block])
    for (g,group) in enumerate(beh)
        # TODO
    end
    replace!(beh.prevPastStopWell, missing=>-1)
    beh = sort(combine(beh, identity),:time)
end

# Checkpoint pre-video data
cp(x, type="arrow") = joinpath(dirname(decode_file), "$(x).$type")
using HDF5
using Arrow
begin
    Arrow.write(cp("cells"), cells)
    Arrow.write(cp("spikes"), spikes)
    Arrow.write(cp("beh"), beh)
    Arrow.write(cp("lfp"), lfp)
    Arrow.write(cp("cycles"), cycles)
    Arrow.write(cp("ripples"), ripples)
    h5open(cp("decode","h5"), "w") do file
        create_group(file, "decode")
        println(file["decode"])
        file["decode/dat"] = dat
        file["decode/x"]   = x
        file["decode/y"]   = y
        file["decode/T"]   = T
    end
    nothing
end


# --------------------
# RUN VIDEO 
# --------------------
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
function select_prob(t, prob::Array{<:Real,3}=dat)
    time  = min(max(t, 1), length(T))
    D = prob[:,:, time]
end
function select_prob4(t, prob::Array{<:Real,4}=dat)
    time  = min(max(t, 1), length(T))
    D = prob[:,:, time, :]
end

## ---------------
## FIGURE WISHLIST
## ---------------
# - Goals
#   - Glow goal
#   - Vector to goal?
# - Sequences
#   - Vector: Start to end
#   - Vector: animal to end
#   - Scatter of maxima
# - Theta
#   - Sweeps suck
# - Ripples
# - Histograms for properties of interest
#   (I.e. best way to test for effects. Example, cosine similarity of animal vec to goal and replay vec)

## FIGURE SETTINGS
visualize = :slider

# Figure
Fig = Figure(resolution=(800,1300))
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
lfp_events = @lift transform($lfp_events, 
                             [:phase,:amp] => ((x,y)->x.*(y.^1.9)) => :phase)

# Grids
# Axes
TT = length(T)
correct = @lift $behavior.correct[1]
@time title_str = @lift "i=$($t), %=$(@sprintf("%2.2f",100*$t/TT)) \ncorr=$($correct), goal=$($behavior.stopWell[1]), vel=$(@sprintf("%2.1f",abs($behavior.velVec[1])))"

axArena  = Axis(Fig[1,1], xlabel="x", ylabel="y", title=title_str)
gNeural = Fig[2,1] = GridLayout(1,1)
axLFPsum = Axis(gNeural[1,1], ylims=(50,150))
axLFP = Axis(gNeural[2:3,1])
axNeural = Axis(gNeural[3:8,1], xlabel="time")
controls = Fig[1:2, 3]


ln_phase_xy = @lift([Point2f(Tuple(x)) for x in eachrow($lfp_events[!,[:time,:phase]])]) 
ln_theta    = lines!(axLFPsum,   ln_phase_xy, color=:white, linestyle=:dash)

# LFP Raw data AXIS
ln_broad_xy = @lift([Point2f(Tuple(x)) for x in eachrow($lfp_events[!,[:time,:broadraw]])]) 
ln_theta_xy = @lift([Point2f(Tuple(x)) for x in eachrow($lfp_events[!,[:time,:raw]])]) 
ln_broad    = lines!(axLFP, ln_broad_xy, color=:gray, linestyle=:dash)
ln_theta    = lines!(axLFP,   ln_theta_xy, color=:white, linestyle=:dash)
axLFP.yticks = 100:100
axLFPsum.yticks = 0:0

# Neural data AXIS
sc_sp_events = @lift([Point2f(Tuple(x)) for x in
                      eachrow($spike_events[!,[:time,:unit]])])
sc = scatter!(axNeural, sc_sp_events, color=:white, markersize=3)
#ln_lfp_phase = @lift([Point2f(Tuple(x)) for x in
#                      eachrow(select($lfp_events, :time, :phase_plot=>x->x.*1))])
#sc = lines!(axNeural, ln_lfp_phase,  color=:gray, linestyle=:dash)
vlines!(axNeural, [0], color=:red, linestyle=:dash)
vlines!(axLFP,    [0], color=:red, linestyle=:dash)
vlines!(axLFPsum, [0], color=:red, linestyle=:dash)

# Arena data AXIS :: Probabilities
#hm_prob = @lift select_prob($t)
#hm = heatmap!(axArena, x,y, hm_prob, colormap=(:bamako,0.8))
theta_prob  = @lift select_prob($t, theta)
ripple_prob = @lift select_prob4($t, ripple)
ca1    = @lift($ripple_prob[:,:,1])
ca1pfc = @lift($ripple_prob[:,:,2])
pfcca1 = @lift($ripple_prob[:,:,3])
pfc    = @lift($ripple_prob[:,:,4])
limits = nanextrema(theta)
hm_theta     = heatmap!(axArena, x,y, theta_prob,   colormap=(:romaO,0.8), colorrange=limits, interpolate=false)
hm_ripple    = heatmap!(axArena, x,y, ca1,  colormap=(:linear_ternary_blue_0_44_c57_n256, 0.8), interpolate=false)
hm_ripple_hp = heatmap!(axArena, x,y, ca1pfc, colormap=(:summer, 0.8), interpolate=false)
hm_ripple_ph = heatmap!(axArena, x,y, pfcca1, colormap=(:spring, 0.8), interpolate=false)
hm_ripple_p  = heatmap!(axArena, x,y, pfc, colormap=(:linear_ternary_red_0_50_c52_n256, 0.8), interpolate=false)
#non_prob = @lift select_prob($t, non)
#hm_non = heatmap!(axArena, x,y, non_prob, colormap=(:bamako,0.8))

# Arena data AXIS :: Behavior
now = 1

point_xy = @lift(Point2f($behavior.x[now], $behavior.y[now]))
line_xy  = @lift([Point2f(b.x, b.y) for b in eachrow($behavior)])
sc_xy    = scatter!(axArena, point_xy, color=:white)
sc_arena = scatter!(axArena, arenaWells.x, arenaWells.y, marker=:star5, markersize=40, color=:gray)
sc_home  = scatter!(axArena, [homeWell.x], [homeWell.y],   marker='𝐇', markersize=25,   color=:gray)
ln_xy = lines!(axArena, line_xy, color=:white, linestyle=:dash)
lines!(axArena, boundary.x, boundary.y, color=:grey)

future₁ = @lift($behavior.stopWell[now])
future₂ = @lift($behavior.futureStopWell[now])
past₁   = @lift($behavior.pastStopWell[now])
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
sc_future₁_well = scatter!(axArena, future₁_well_xy, alpha=0.5, marker='𝐅', markersize=60, color=correctColor, glowwidth=29)
sc_future₂_well = scatter!(axArena, future₂_well_xy, alpha=0.5, marker='𝐟', markersize=60, color=:white, glowwidth=5)
sc_past_well    = scatter!(axArena, past₁_well_xy,    alpha=0.5, marker='𝐏', markersize=60, color=:white, glowwidth=5)
if doPrevPast
    past₂   = @lift($behavior.pastStopWell[now])
    past₂_well_xy    = @lift (($past₂ == -1) || ($past₂ == $future₁)) ? Point2f(NaN, NaN) : Point2f(wells.x[$past₂], wells.y[$past₂])
end

function play_graphic(stop=length(T))
    framerate = 90
    start = t[]
    timestamps = range(start, stop, step=1)
    P=Progress(stop-start, desc="Video")
    #recording = plotsdir("mpp_decode", "withBehVideo=$usevideo", outputVideo)
    path = joinpath(dirname(decode_file), "decode_split=$(split_num)_start=$(start)_stop=$stop.mp4")
    record(Fig, path, timestamps; framerate=framerate, compression=10) do stamp
        t[] = stamp
        next!(P)
    end
end


# -------------
# VISUALIZE
# -------------
if visualize == :video
    play_graphic()
elseif visualize == :slider
end

