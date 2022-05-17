ripples, cycles = annotate_vector_info(ripples, cycles, beh, lfp, dat, x, y, T)

# Cast to Float32(for makie) and range for plotting
lfp.phase_plot = utils.norm_extrema(lfp.phase, extrema(spikes.unit))
@assert lfp != nothing
lfp.phase = utils.norm_extrema(lfp.phase, (-pi,pi))
@assert lfp != nothing
lfp.raw      = Float32.(utils.norm_extrema(lfp.raw,      extrema(spikes.unit)))
@assert lfp != nothing
lfp.broadraw = Float32.(utils.norm_extrema(lfp.broadraw, extrema(spikes.unit)))
@assert lfp != nothing
lfp.cycle, lfp.phase = Int32.(lfp.cycle), Float32.(lfp.phase)
@assert lfp != nothing
cycles = cycles[cycles.cycle.!=-1,:]


