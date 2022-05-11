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

