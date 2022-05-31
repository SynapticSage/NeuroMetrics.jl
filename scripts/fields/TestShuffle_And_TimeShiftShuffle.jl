quickactivate("/home/ryoung/Projects/goal-code/")
@time include(scriptsdir("fields", "Initialize.jl"))
@time include(scriptsdir("fields", "TimeShift_setsOfInterest.jl"))
includet(srcdir("shuffle.jl"))
info      = field.info
timeshift = field.timeshift

_, spikes = raw.register(beh, spikes; transfer=["velVec"], on="time")

import Base.Threads: @spawn
using ThreadSafeDicts
using .timeshift
sp = copy(spikes)
convertToMinutes, runshifts = true, false

# -----------------------
# Testing shuffle methods
# -----------------------
# @time shuffle.jitterBy(spikes, distribution=:uniform, width=:traj,    data=beh) # 2.2 seconds
# @time shuffle.jitterBy(spikes, distribution=:uniform, width=:session, data=beh) # 5 seconds
# 
# @time shuffle.jitterAllSpikes(spikes; distribution=:uniform, width=:traj, data=beh) # 0.36 second
# @time shuffle.jitterAllSpikes(spikes; distribution=:uniform, width=:session, data=beh) # 0.13 seconds
# 
# -----------------------------
# Testing time-shifted shuffles
# -----------------------------
props = ["x", "y"]
splitby=["unit", "area"]
kws=(;splitby, filters=merge(filt.speed_lib, filt.cellcount))
newkws = (; kws..., resolution=40, gaussian=2.3*0.5, props=props,
          filters=merge(kws.filters))
shifts = -4:0.2:4

if convertToMinutes
    @info "Coverting shifts to minutes"
    shifts = shifts./60
    #using Serialization
    #safe_dict = deserialize(datadir("safe_dict.serial"))
else
    safe_dict = ThreadSafeDict()
end

result = @time timeshift.get_field_shift(beh, spikes, shifts;
                         multi=:single, postfunc=info.information, 
                         newkws...) # 2 minutes, roughly

@assert sp.time == spikes.time

shuf_result = @time timeshift.get_field_shift_shuffles(beh, spikes, shifts; # 4-16 hours
                         multi=:single, postfunc=info.information, 
                         shuffle_func=shuffle.by,
                         exfiltrateAfter=Inf,
                         safe_dict=safe_dict,
                         newkws...)

timeshift.saveshifts(result, shuf_result, shifts=shifts, metric="info", 
                     fieldkws=newkws)

out, cv = timeshift.get_field_crossval_sk(beh, spikes, shifts; n_folds=2, 
                             multi=:single, postfunc=field.info.information, 
                             newkws...)

# TESTING CORRECTION METHODS
D = timeshift.loadshifts(;shifts=shifts, metric="info", newkws...)
main, shuf= info_to_dataframe(D[:main]), info_to_dataframe(D[:shuffle])

