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

load_from_checkpoint          = false
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

@showprogress 0.1 "Processing splits" for (split_num, split_type) in collect(Iterators.product([0,1,2,3], ["test",]))

    video="/Volumes/Colliculus/RY16_experiment/actualVideos/RY16_69_0$(epoch)_CMt.1.mp4"
    outputVideo = "animation.$(decoder_type)_$(transition_type)_$(split_type)_$(split_num)_$(variable)_$(basename(video))"
    decode_file = raw.decodepath(animal, day, epoch, transition="empirical",
                                 method="sortedspike", split=split_num,
                                 type=split_type, speedup=20.0)
    if usevideo
        stream = VideoIO.open(video)
        vid = VideoIO.openvideo(stream)
    end

    @infiltrate
    h5file = decode.pathname("split=$(split_num)_decode", decode_file, "h5")
    if isfile(h5file)
        @infiltrate
        @info "$h5file is already a file...skipping"
        continue
    end

    # -----------
    # DATASETS
    # -----------
    global beh, cells, spikes, lfp, ripples, cycles, wells, x, y, T, dat
    lfp = beh = ripples = spikes = dat = non = ripple = theta = nothing
    beh    = raw.load_behavior(animal, day)
    spikes = raw.load_spikes(animal,   day)
    cells  = raw.load_cells(animal,    day)
    lfp    = @subset(raw.load_lfp(animal, day),
                           :tetrode.==ca1_tetrode)[!,[:time, :raw, :phase,
                                                      :amp, :broadraw]]
    task   = raw.load_task(animal,     day)
    ripples= raw.load_ripples(animal,  day)
    wells = task[(task.name.=="welllocs") .& (task.epoch .== epoch), :]
    D = raw.load_decode(decode_file)
    @info "Loaded all data"

    # Normalize times
    mint = minimum(beh[beh.epoch.==epoch,:].time)
    beh, spikes, lfp, ripples, D = raw.normalize_time(beh, spikes, lfp, ripples, D; 
                                                      timefields=Dict(4=>["start", "stop", "time"]));
    x, y, T, dat = D["x_position"], D["y_position"], D["time"], D[variable]
    dat = permutedims(dat, [2,1,3])
    dat = decode.quantile_threshold(dat, thresh[variable])
    D = nothing
    nmint = minimum(beh[beh.epoch.==epoch,:].time)
    @info "Initial munge"

    if remove_nonoverlap
        spikes, beh, lfp, ripples, T_inds = 
             raw.keep_overlapping_times(spikes, beh, lfp, ripples, T;
                                        returninds=[5])
        dat, T = dat[:,:,T_inds], T[T_inds]
        @info "Removed overlap"
    end
    [extrema(x.time) for x in (lfp, spikes, ripples, beh)]

    # --------------------
    # Preprocess BEHAHVIOR
    # TODO Turn this into a function
    # --------------------
    beh = annotate_pastFutureGoals(beh; doPrevPast=false)
    @info "Added behavioral fields"

    # --------------
    # Preprocess LFP
    # --------------
    ripples      = velocity_filter_ripples(beh,ripples)
    lfp, cycles  = get_theta_cycles(lfp, beh)
    lfp          = curate_lfp_theta_cycle_and_phase(lfp, cycles)
    lfp, ripples = annotate_ripples_to_lfp(lfp, ripples)
    ripples, cycles = annotate_vector_info(ripples, cycles, beh, lfp, dat, x, y, T)

    # Cast to Float32(for makie) and range for plotting
    lfp.phase_plot = utils.norm_extrema(lfp.phase, extrema(spikes.unit))
    lfp.phase = utils.norm_extrema(lfp.phase, (-pi,pi))
    lfp.raw      = Float32.(utils.norm_extrema(lfp.raw,      extrema(spikes.unit)))
    lfp.broadraw = Float32.(utils.norm_extrema(lfp.broadraw, extrema(spikes.unit)))
    lfp.cycle, lfp.phase = Int32.(lfp.cycle), Float32.(lfp.phase)
    cycles = cycles[cycles.cycle.!=-1,:]
    @info "Added lfp, cycle, and ripple fields"

    global theta, ripple, non
    lfp, theta, ripples, non = begin
        # ------------------------------
        # Separate DECODES by EVENTS
        # ------------------------------
        theta, ripple, non = separate_theta_ripple_and_non_decodes(T, lfp, dat;
                                                                    doRipplePhase=doRipplePhase)
        if dosweep
            theta, ripples = convert_to_sweeps(lfp, theta, ripples;
                                               doRipplePhase=doRipplePhase)
        end
        no_cycle_happening = utils.squeeze(all(isnan.(theta), dims=(1,2)))
        lfp[!, :color_phase] = RGBA.(get(ColorSchemes.romaO, lfp[!, :phase]))
        lfp[utils.searchsortednearest.([lfp.time], T[no_cycle_happening]), :color_phase] .= RGBA{Float64}(0, 0, 0, 0)
        lfp, theta, ripples, non
    end
    @info "Split decode into lfp, theta, and non"

    # Checkpoint pre-video data
    @info "Saving"
    decode.save_checkpoint(Main, decode_file; split=split_num, overwrite=false)
    exit()
end

