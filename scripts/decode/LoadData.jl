import NetCDF
import Load


# -----------
# DATASETS
# -----------
global beh, cells, spikes, lfp, ripples, cycles, wells, x, y, T, dat
lfp = beh = ripples = spikes = dat = non = ripple = theta = nothing
GC.gc()

#pfc_tetrodes = reshape([36, 39, 40, 47, 48], 1, 5)
pfc_tetrodes = []

@time beh    = Load.load_behavior(animal, day)
@time spikes = Load.load_spikes(animal,   day)
@time cells  = Load.load_cells(animal,    day)
@time task   = Load.load_task(animal,     day)
@time ripples= Load.load_ripples(animal,  day)
@time D      = Load.load_decode(decode_file)
@time lfp    = Load.load_lfp(animal, day;
							tet=[ca1_tetrode, pfc_tetrodes...],
							vars=[:time, :raw, :phase, :amp, :broadraw, :tetrode])
#beh, spikes, cells ,task, ripples, D = fetch(beh), fetch(spikes), fetch(cells), 
#                                fetch(task), fetch(ripples), fetch(D)


if !isempty(pfc_tetrodes)
    pfc_lfp = @subset(lfp, utils.squeeze(any(:tetrode .== pfc_tetrodes, dims=2)))
    lfp = @subset(lfp, :tetrode.==ca1_tetrode)
    GC.gc()
end

@info "Loaded all data"

wells = task[(task.name.=="welllocs") .& (task.epoch .== epoch), :]
homeWell, arenaWells = begin
    hw = argmax(StatsBase.fit(StatsBase.Histogram, filter(b->b!=-1,beh.stopWell), 1:6).weights)
    hw=5
    @info "Homewell = $hw"
    wells[hw,:], wells[setdiff(1:5,hw),:]
end
boundary = task[(task.name.=="boundary") .& (task.epoch .== epoch), :]
append!(boundary, DataFrame(boundary[1,:]))

# Normalize times
if lfp === nothing
	@error "fuck"
end

mint = minimum(beh[beh.epoch.==epoch,:].time)
if !isempty(pfc_tetrodes)
    beh, spikes, lfp, pfc_lfp, ripples, D = Load.normalize_time(beh, spikes, lfp, pfc_lfp, ripples, D; 
                                                      timefields=Dict(5=>["start", "stop", "time"]));
else
    beh, spikes, lfp, ripples, D = Load.normalize_time(beh, spikes, lfp, ripples, D; 
                                                      timefields=Dict(4=>["start", "stop", "time"]));
end
x, y, T, dat = D["x_position"], D["y_position"], D["time"], D[variable]
dat = permutedims(dat, [2,1,3])

if dothresh
    @info "Appliying thresh=$thresh in LoadData.jl"
    dat = decode.quantile_threshold(dat, thresh[variable])
end
D = nothing

nmint = minimum(beh[beh.epoch.==epoch,:].time)

@info "Initial lfp cycle steps"
ripples      = velocity_filter_ripples(beh,ripples)
@assert ripples isa DataFrame
lfp, cycles  = get_theta_cycles(lfp, beh)
lfp          = curate_lfp_theta_cycle_and_phase(lfp, cycles)
lfp, ripples = annotate_ripples_to_lfp(lfp, ripples)

if lfp === nothing
	@error "fuck"
end
if cycles === nothing
	@error "fuck"
end


ripples_copy = copy(ripples)
if remove_nonoverlap
    spikes, beh, lfp, pfc_lfp, ripples, T_inds, cycles = 
         raw.keep_overlapping_times(spikes, beh, lfp, pfc_lfp, ripples, T, cycles;
                                    returninds=[6])
    dat, T = dat[:,:,T_inds], T[T_inds]
    @info "Removed overlap"
end
[extrema(x.time) for x in (lfp, spikes, ripples, beh)]
if lfp === nothing
	@error "fuck"
end
if cycles === nothing
	@error "fuck"
end
@assert !(isempty(beh))
