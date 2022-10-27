quickactivate("/home/ryoung/Projects/goal-code/")
@time include(scriptsdir("fields", "Initialize.jl"))
@time include(scriptsdir("fields", "TimeShift_setsOfInterest.jl"))
includet(srcdir("shuffle.jl"))
info      = field.info
timeshift = field.timeshift

nbins = 50
raw.behavior.annotate_relative_xtime!(beh);
beh.trajreltime_bin = floor.(beh.trajreltime * (nbins-1));
_, spikes = raw.register(beh, spikes; 
                         transfer=["trajreltime","trajreltime_bin","velVec"], 
                         on="time");

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
#
spikesnew = @time shuffle.permuteBy(spikes; split=:trajreltime_bin, sort=false, keepold=true)
for f in ("trajreltime","trajreltime_bin")
    spikesnew[:,f*"old"] = spikesnew[:,f]
end
_, spikesnew = raw.register(beh, spikesnew; 
                         transfer=["trajreltime","trajreltime_bin","velVec"], 
                         on="time");
#sort!(spikesnew, :time)
@df spikesnew[1:1000:end, :] plot(scatter(:timeold, :time, alpha=0.2, xlabel="old",ylabel="new"), histogram(:time .- :timeold, xlabel="difference"))
@df spikesnew[1:1000:end,:]  plot(
                                 scatter(:trajreltime, :trajreltimeold, title="Tiny bit of jitter\nunderlying reltime", xlabel="new",ylabel="old",alpha=0.5), 
                                 scatter(:trajreltime_bin, :trajreltime_binold, title="Preserved bin", xlabel="new",ylabel="old",alpha=0.5))

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
                         preset=:cDt_t,
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

