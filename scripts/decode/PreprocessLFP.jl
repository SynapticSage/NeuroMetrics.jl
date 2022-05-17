ripples, cycles = annotate_vector_info(ripples, cycles, beh, lfp, dat, x, y, T)

# Cast to Float32(for makie) and range for plotting
lfp.phase_plot = utils.norm_extrema(lfp.phase, extrema(spikes.unit))
lfp.phase = utils.norm_extrema(lfp.phase, (-pi,pi))
lfp.raw      = Float32.(utils.norm_extrema(lfp.raw,      extrema(spikes.unit)))
lfp.broadraw = Float32.(utils.norm_extrema(lfp.broadraw, extrema(spikes.unit)))
lfp.cycle, lfp.phase = Int32.(lfp.cycle), Float32.(lfp.phase)
cycles = cycles[cycles.cycle.!=-1,:]


