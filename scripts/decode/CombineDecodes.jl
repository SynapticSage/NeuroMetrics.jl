# -------
using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
include(scriptsdir("scripts","decode","Initialize.jl"))

@showprogress 0.1 "Processing splits" for (split_num, split_type) in collect(Iterators.product([0,1,2,3], ["test",]))

    decode_file = raw.decodepath(animal, day, epoch, transition="empirical",
                                 method="sortedspike", split=split_num,
                                 type=split_type, speedup=20.0)
    if usevideo
        stream = VideoIO.open(video)
        vid = VideoIO.openvideo(stream)
    end

    h5file = decode.pathname("split=$(split_num)_decode", decode_file, "h5")
    if isfile(h5file) && split_num != 0
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
    decode.save_checkpoint(Main, decode_file; split=split_num, overwrite=true)
    break
end

