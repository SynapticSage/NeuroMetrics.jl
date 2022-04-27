quickactivate("/home/ryoung/Projects/goal-code/")
@time include(scriptsdir("fields", "Initialize.jl"))
includet(srcdir("shuffle.jl"))
includet(srcdir("field/info.jl"))
includet(srcdir("field/timeshift.jl"))
_, spikes = raw.register(beh, spikes; transfer=["velVec"], on="time")
import Base.Threads: @spawn
sp = copy(spikes)

# Testing shuffle methods
# -----------------------
@time shuffle.by(spikes, distribution=:uniform, width=:traj,    data=beh) # 2.2 seconds
@time shuffle.by(spikes, distribution=:uniform, width=:session, data=beh) # 2.9 seconds

@time shuffle.byspike(spikes; distribution=:uniform, width=:traj, data=beh) # 0.36 seconds
@time shuffle.byspike(spikes; distribution=:uniform, width=:session, data=beh) # 0.13 seconds

# Testing time-shifted shuffles
# -----------------------------
props = ["x", "y"]
splitby=["unit", "area"]
kws=(;splitby, filters=merge(filt.speed_lib, filt.cellcount))
newkws = (; kws..., resolution=40, gaussian=2.3*0.5, props=props,
          filters=merge(kws.filters))
shifts = -4:0.2:4

# Single thread
#
result = @time timeshift.get_field_shift(beh, spikes, shifts;
                         multi=:single, postfunc=info.information, 
                         newkws...) # 2 minutes, roughly

@assert sp.time == spikes.time

shuf_result = @time timeshift.get_field_shift_shuffles(beh, spikes, shifts; # 4-16 hours
                         multi=:single, postfunc=info.information, 
                         shuffle_func=shuffle.by,
                         exfiltrateAfter=25,
                         newkws...)

# Distributed
using Distributed
using Dagger
addprocs(8, enable_threaded_blas=true)
@everywhere include("/home/ryoung/Projects/goal-code/scripts/fields/Include.jl")
@time shuf_result = timeshift.get_field_shift_shuffles(beh, spikes, shifts; # 16 hours
                         multi=:distributed, postfunc=info.information, 
                         shuffle_func=shuffle.by,
                         newkws...)

# Threading
Threads.nthreads() = 8
@time shuf_result = timeshift.get_field_shift_shuffles(beh, spikes, shifts; # 16 hours
                         multithread=true, postfunc=info.information, 
                         shuffle_func=shuffle.by,
                         exfiltrateAfter=25,
                         newkws...)
