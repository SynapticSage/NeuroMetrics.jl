# -------
# IMPORTS
# -------
using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
using Pushover, Revise, Interact, Blink, Mux, ProgressMeter
using Statistics, NaNStatistics
using VideoIO
using GLMakie
using ColorSchemes, Colors
import ColorSchemeTools 
using DataFrames, DataFramesMeta
using Printf
using StatsPlots: @df
import StatsBase
using Infiltrator
set_theme!(theme_dark())
__revise_mode__ = :eval
includet(srcdir("raw.jl"))
includet(srcdir("table.jl"))
includet(srcdir("utils.jl"))
includet(srcdir("decode.jl"))
import .raw, .table, .utils
using .decode

# Debugger?
ENV["JULIA_DEBUG"] = nothing

# -----------------
# HELPER FUNCTIONS 
# -----------------
na = [CartesianIndex()]
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

# -----------
# Decode file locations
# -----------
dir = plotsdir("mpp_decode", "withBehVideo=$usevideo")
if !(isdir(dir))
    mkdir(dir)
end
(split_num, split_type) = collect(Iterators.product([0,1,2], 
                                                ["test","train"]))[1]
video="/Volumes/Colliculus/RY16_experiment/actualVideos/RY16_69_0$(epoch)_CMt.1.mp4"
outputVideo = "animation.$(decoder_type)_$(transition_type)_$(split_type)_$(split_num)_$(variable)_$(basename(video))"
decode_file = raw.decodepath(animal, day, epoch, transition="empirical",
                             method="sortedspike", split=split_num,
                             type=split_type, speedup=20.0)
stream = VideoIO.open(video)
vid = VideoIO.openvideo(stream)

if !(load_from_checkpoint)
    # -----------
    # DATASETS
    # -----------
    @time beh    = raw.load_behavior(animal, day)
    @time spikes = raw.load_spikes(animal,   day)
    @time cells  = raw.load_cells(animal,    day)
    @time lfp    = @subset(raw.load_lfp(animal, day),
                           :tetrode.==ca1_tetrode)[!,[:time, :raw, :phase, :amp, :broadraw]]
    @time task   = raw.load_task(animal,     day)
    @time ripples= raw.load_ripples(animal,  day)
    wells = task[(task.name.=="welllocs") .& (task.epoch .== epoch), :]
    D = raw.load_decode(decode_file)

    # Normalize times
    mint = minimum(beh[beh.epoch.==epoch,:].time)
    beh, spikes, lfp, ripples, D = raw.normalize_time(beh, spikes, lfp, ripples, D; 
                                                      timefields=Dict(4=>["start", "stop", "time"]));
    x, y, T, dat = D["x_position"], D["y_position"], D["time"], D[variable]
    dat = permutedims(dat, [2,1,3])
    dat = decode.quantile_threshold(dat, thresh[variable])
    D = nothing
    nmint = minimum(beh[beh.epoch.==epoch,:].time)

    if remove_nonoverlap
        spikes, beh, lfp, ripples, T_inds = raw.keep_overlapping_times(spikes, beh, lfp, ripples, T;
                                                                      returninds=[5])
        dat, T = dat[:,:,T_inds], T[T_inds]
    end
    [extrema(x.time) for x in (lfp, spikes, ripples, beh)]

    # --------------------
    # Preprocess BEHAHVIOR
    # TODO Turn this into a function
    # --------------------
    beh = annotate_pastFutureGoals(beh)

    # --------------
    # Preprocess LFP
    # --------------
    ripples      = velocity_filter_ripples(beh,ripples)
    lfp, cycles  = get_theta_cycles(lfp)
    lfp          = curate_lfp_theta_cycle_and_phase(lfp, cycles)
    lfp, ripples = annotate_ripples_to_lfp(lfp, ripples)
    ripples, cycles = decode.annotate_vector_info(ripples, cycles, beh, dat, x, y, T)

    # Cast to Float32(for makie) and range for plotting
    lfp.phase_plot = utils.norm_extrema(lfp.phase, extrema(spikes.unit))
    lfp.phase = utils.norm_extrema(lfp.phase, (-pi,pi))
    lfp.raw      = Float32.(utils.norm_extrema(lfp.raw,      extrema(spikes.unit)))
    lfp.broadraw = Float32.(utils.norm_extrema(lfp.broadraw, extrema(spikes.unit)))
    lfp.cycle, lfp.phase = Int32.(lfp.cycle), Float32.(lfp.phase)

    # Checkpoint pre-video data
    save_checkpoint(Main)
else
    D = decode.load_checkpoint(Main, decode_file)
    for (key,value) in D
        eval(Meta.parse("$key = D[:$key]"))
    end
    fix_complex(x) = x.re + (x.im)im
    beh.velVec = fix_complex.(beh.velVec)
    ripples.velVec = fix_complex.(ripples.velVec)
    lfp.velVec = fix_complex.(lfp.velVec)
end

cycles = cycles[cycles.cycle.!=-1,:]
lfp, theta, ripples, non = begin
    # ------------------------------
    # Separate DECODES by EVENTS
    # ------------------------------
    theta, ripple, non = separate_theta_ripple_and_non_decodes(T, lfp, dat;
                                                                doRipplePhase=doRipplePhase)
    frac_print(frac) = @sprintf("Fraction times = %2.4f",frac)
    frac = mean(all(isnan.(theta) .&& isnan.(ripple), dims=(1,2)))
    frac_print(frac)
    frac_rip = vec(mean(all(isnan.(ripple),dims=(1,2)), dims=(1,2,3)))
    frac_print.(frac_rip)
    if dosweep
        theta, ripples = convert_to_sweeps(lfp, theta, ripples;
                                           doRipplePhase=doRipplePhase)
    end
    no_cycle_happening = utils.squeeze(all(isnan.(theta), dims=(1,2)))
    lfp[!, :color_phase] = RGBA.(get(ColorSchemes.romaO, lfp[!, :phase]))
    lfp[utils.searchsortednearest.([lfp.time], T[no_cycle_happening]), :color_phase] .= RGBA{Float64}(0, 0, 0, 0)
    lfp, theta, ripples, non
end

homeWell, arenaWells = begin
    hw = argmax(StatsBase.fit(StatsBase.Histogram, filter(b->b!=-1,beh.stopWell), 1:6).weights)
    @info "Homewell = $hw"
    wells[hw,:], wells[setdiff(1:5,hw),:]
end
boundary = task[(task.name.=="boundary") .& (task.epoch .== epoch), :]
append!(boundary, DataFrame(boundary[1,:]))

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

## ---------------
## FIGURE WISHLIST
## ---------------
# - Goals
#   - Glow goal
#   - □ Vector to goal?
# - Errors
#   - ✔ Show actual poked well that caused the error
#   - □ Pass poke data over trial in seperate column
# - Sequences
#   - ✓ Vector: Start to end
#   - □ Vector: animal to end
#   - □ Scatter of maxima
# - Theta
#   □ Sweeps suck
#   □ θ sweep arrow : start at trough, end at peak
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
@time behavior     = @lift select_range($t, T, data=beh, Δ_bounds=[-0.01, 0.60])
@time spike_events = @lift select_range($t, T, data=spikes, Δ_bounds=Δ_bounds)
@time lfp_events   = @lift select_est_range($t, T, tr["lfp"], Δt["lfp"], 
                                            data=lfp, Δ_bounds=Δ_bounds)
@time ripple_events = @lift decode.select_events($t, T, events=ripples)
@time cycle_events  = @lift decode.select_events($t, T, events=cycles)
theta_prob  = @lift select_prob($t, T; prob=theta)
ripple_prob = @lift select_prob4($t, T; prob=ripple)

lfp_now    = @lift Int32(round(size($lfp_events,1)/2))
lfp_events = @lift transform($lfp_events, 
                             [:phase,:amp] => ((x,y)->x.*(y.^1.9)) => :phaseamp)

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

ln_theta_color = @lift mean($lfp_events.color_phase[-10+$lfp_now:10+$lfp_now])
ln_phase_xy = @lift([Point2f(Tuple(x)) for x in eachrow($lfp_events[!,[:time,:phaseamp]])]) 
ln_theta    = lines!(axLFPsum,   ln_phase_xy, color=ln_theta_color, linestyle=:dash)

# LFP Raw data AXIS
ln_broad_xy = @lift([Point2f(Tuple(x)) for x in eachrow($lfp_events[!,[:time,:broadraw]])]) 
ln_theta_xy = @lift([Point2f(Tuple(x)) for x in eachrow($lfp_events[!,[:time,:raw]])]) 
ln_broad    = lines!(axLFP, ln_broad_xy, color=:gray,  linestyle=:dash)
ln_theta    = lines!(axLFP, ln_theta_xy, color=ln_theta_color, linestyle=:dash)
axLFP.yticks = 100:100
axLFPsum.yticks = 0:0

# Neural data AXIS
colorby = nothing
spike_colors = @lift begin
    if colorby == nothing
        :white
    else
        $behavior[!, colorby]
    end
end
spike_colorrange = begin
    if colorby == nothing
        (-Inf, Inf)
    else
        extrema(beh[!, colorby])
    end
end
sc_sp_events = @lift([Point2f(Tuple(x)) for x in
                      eachrow($spike_events[!,[:time,:unit]])])
sc = scatter!(axNeural, sc_sp_events, color=spike_colors, markersize=3, colorrange=spike_colorrange)
#ln_lfp_phase = @lift([Point2f(Tuple(x)) for x in
#                      eachrow(select($lfp_events, :time, :phase_plot=>x->x.*1))])
#sc = lines!(axNeural, ln_lfp_phase,  color=:gray, linestyle=:dash)
vlines!(axNeural, [0], color=:red, linestyle=:dash)
vlines!(axLFP,    [0], color=:red, linestyle=:dash)
vlines!(axLFPsum, [0], color=:red, linestyle=:dash)

# Arena data AXIS :: Probabilities
#hm_prob = @lift select_prob($t)
#hm = heatmap!(axArena, x,y, hm_prob, colormap=(:bamako,0.8))
ca1    = @lift($ripple_prob[:,:,1])
ca1pfc = @lift($ripple_prob[:,:,2])
pfcca1 = @lift($ripple_prob[:,:,3])
pfc    = @lift($ripple_prob[:,:,4])
limits = nanextrema(theta)
hm_theta     = heatmap!(axArena, x,y, theta_prob, colormap=(:romaO,0.8), colorrange=limits, interpolate=false)
hm_ripple    = heatmap!(axArena, x,y, ca1,        colormap=(:linear_ternary_blue_0_44_c57_n256, 0.8), interpolate=false)
hm_ripple_hp = heatmap!(axArena, x,y, ca1pfc,     colormap=(:summer, 0.8), interpolate=false)
hm_ripple_ph = heatmap!(axArena, x,y, pfcca1,     colormap=(:spring, 0.8), interpolate=false)
hm_ripple_p  = heatmap!(axArena, x,y, pfc,        colormap=(:linear_ternary_red_0_50_c52_n256, 0.8), interpolate=false)
if plotNon
    non_prob = @lift select_prob($t, T, prob=non)
    hm_non = heatmap!(axArena, x,y, non_prob, colormap=(:bamako,0.8))
end

# Arena data AXIS :: Behavior
now = 1

# ----------------
# Animal and wells
# ----------------
point_xy = @lift(Point2f($behavior.x[now], $behavior.y[now]))
line_xy  = @lift([Point2f(b.x, b.y) for b in eachrow($behavior)])
sc_xy    = scatter!(axArena, point_xy, color=:white)
sc_arena = scatter!(axArena, arenaWells.x, arenaWells.y, marker=:star5, markersize=40, color=:gray)
sc_home  = scatter!(axArena, [homeWell.x], [homeWell.y],   marker='𝐇', markersize=25,   color=:gray)
ln_xy = lines!(axArena, line_xy, color=:white, linestyle=:dash)
lines!(axArena, boundary.x, boundary.y, color=:grey)

# ----------------
# Sequence vectors
# ----------------
cycle_start_xy  = @lift isempty($cycle_events) ? [NaN, NaN]  : $cycle_events.dec₀
cycle_uv        = @lift isempty($cycle_events) ? [NaN, NaN]  : $cycle_events.dec₀₁
cycle_arrows = arrows!(axArena, cycle_start_xy[1], cycle_start_xy[2],
                                cycle_uv[1], cycle_uv[2])

ripple_start_xy  = @lift isempty($ripple_events) ? [NaN, NaN]  : $ripple_events.dec₀
ripple_uv        = @lift isempty($ripple_events) ? [NaN, NaN]  : $ripple_events.dec₀₁
ripple_arrows = arrows!(axArena, ripple_start_xy[1], ripple_start_xy[2],
                                 ripple_uv[1], ripple_uv[2])

# ---------------
# Sequence maxima
# ---------------

# ----------------
# Goal vectors
# ----------------

# ---------------
# Goals and pokes
# ---------------
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
sc_future₁_well = scatter!(axArena, future₁_well_xy,  alpha=0.5, marker='𝐅', markersize=60, color=correctColor, glowwidth=29)
sc_future₂_well = scatter!(axArena, future₂_well_xy,  alpha=0.5, marker='𝐟', markersize=60, color=:white,       glowwidth=5)
sc_past_well    = scatter!(axArena, past₁_well_xy,    alpha=0.5, marker='𝐏', markersize=60, color=:white,       glowwidth=5)
if doPrevPast
    past₂   = @lift($behavior.pastStopWell[now])
    past₂_well_xy    = @lift (($past₂ == -1) || ($past₂ == $future₁)) ? Point2f(NaN, NaN) : Point2f(wells.x[$past₂], wells.y[$past₂])
end

# Pokes if error
poke_1_xy = @lift (!($correct==0) && $behavior.poke_1[now]==1) ? Point2f(wells.x[1], wells.y[1]) : Point2f(NaN, NaN)
poke_2_xy = @lift (!($correct==0) && $behavior.poke_2[now]==1) ? Point2f(wells.x[2], wells.y[2]) : Point2f(NaN, NaN)
poke_3_xy = @lift (!($correct==0) && $behavior.poke_3[now]==1) ? Point2f(wells.x[3], wells.y[3]) : Point2f(NaN, NaN)
poke_4_xy = @lift (!($correct==0) && $behavior.poke_4[now]==1) ? Point2f(wells.x[4], wells.y[4]) : Point2f(NaN, NaN)
poke_5_xy = @lift (!($correct==0) && $behavior.poke_5[now]==1) ? Point2f(wells.x[5], wells.y[5]) : Point2f(NaN, NaN)
poke_1_sc = scatter!(axArena, poke_1_xy, marker='↑', markersize=30, color=:white, glow_width=10)
poke_2_sc = scatter!(axArena, poke_2_xy, marker='↑', markersize=30, color=:white, glow_width=10)
poke_3_sc = scatter!(axArena, poke_3_xy, marker='↓', markersize=30, color=:white, glow_width=10)
poke_4_sc = scatter!(axArena, poke_4_xy, marker='↑', markersize=30, color=:white, glow_width=10)
poke_5_sc = scatter!(axArena, poke_5_xy, marker='←', markersize=30, color=:white, glow_width=10)

function play_graphic(stop=length(T)-2000; framerate = 90)
    start = t[]
    timestamps = range(start, stop, step=1)
    P=Progress(stop-start, desc="Video")
    path = joinpath(dirname(decode_file), "decode_split=$(split_num)_start=$(start)_stop=$stop.mp4")
    record(Fig, path, timestamps; framerate=framerate, compression=10) do stamp
        try
            t[] = stamp
        catch
            @warn "timestamp failed with t=$(t[])"
        end
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

