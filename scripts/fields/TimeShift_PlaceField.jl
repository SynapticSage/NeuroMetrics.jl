quickactivate("/home/ryoung/Projects/goal-code/")
using Base.Threads: @spawn
using DataFrames
using StatsPlots
using Statistics
using StatsPlots
using DataFrames
includet(srcdir("field.jl"))
includet(srcdir("filt.jl"))
includet(srcdir("raw.jl"))
includet(srcdir("field/timeshift.jl"))
includet(srcdir("field/info.jl"))
includet(srcdir("utils.jl"))
includet(srcdir("table.jl"))
@time spikes, beh = raw.load("RY16", 36, data_source=["spikes","behavior"])

correctToMinutes = true
if correctToMinutes
    c(x) = x./60
else
    c(x) = identity(x)
end

# MULTIPLE TIME SHIFT NOTES
#
# -----------
# 16 threads
# -----------
# (run 1) 400 seconds for 20, with multi-threading, with 90% of that compile time
# (run 2) 180 seconds 85% compile time
# (run 3) 131 seconds 76% compile time
# (run 4) 131 seconds 76% compile time
#
# ---------
# 4 threads
# ---------
#
# ---------
# 1 thread
# ---------
# (run 1) 108 seconds for 20 shifts without multi-threading, with 10% compile time
#
# ------------
# More shifts
# ------------
# For 200 shifts, we're in the domain of 30 minutes
#
#
# A (( ---- PLACE ---- ))
props = ["x", "y"]
splitby=["unit", "area"]
kws=(;splitby, filters=merge(filt.speed_lib, filt.cellcount), multi=:single, postfunc=info.information)
newkws = (; kws..., resolution=40, gaussian=2.3*0.5, props=props, filters=merge(kws.filters))

place = @spawn @time timeshift.get_field_shift(beh, spikes, c(-2:0.1:2); newkws...);
place_fineNarrow = @spawn @time timeshift.get_field_shift(beh, spikes, c(-2:0.01:2); newkws...);
place_broadTimes = @spawn @time timeshift.get_field_shift(beh, spikes, c(-4:0.1:4); newkws...);

newkws = (; kws..., resolution=40, gaussian=2.3*0.5, props=props, filters=merge(kws.filters, filt.correct))
place_correctOnly = @time timeshift.get_field_shift(beh, spikes, c(-2:0.1:2); newkws...);

newkws = (; kws..., resolution=40, gaussian=2.3*0.5, props=props, filters=merge(kws.filters, filt.incorrect))
place_incorrectOnly = @spawn @time timeshift.get_field_shift(beh, spikes, c(-2:0.1:2); newkws...);

newkws = (; kws..., resolution=40, gaussian=2.3*0.5, props=props, filters=merge(kws.filters, filt.nontask))
place_nontaskOnly = @spawn @time timeshift.get_field_shift(beh, spikes, c(-2:0.1:2); newkws...);

newkws = (; kws..., resolution=40, gaussian=2.3*0.5, props=props, filters=merge(kws.filters, filt.cue, filt.correct))
place_correctCue = @spawn @time timeshift.get_field_shift(beh, spikes, c(-2:0.1:2); newkws...);

newkws = (; kws..., resolution=40, gaussian=2.3*0.5, props=props, filters=merge(kws.filters, filt.mem, filt.correct))
place_correctMem = @spawn @time timeshift.get_field_shift(beh, spikes, c(-2:0.1:2); newkws...);


# B (( ---- GOAL ---- ))
filters = merge(kws.filters,
                filt.correct,
                filt.notnan("currentAngle"), 
                filt.minmax("currentPathLength", 2, 150))
props = ["currentAngle", "currentPathLength"]
newkws = (; kws..., filters=merge(kws.filters, filters),
          resolution=40, gaussian=2.3*0.5, props=props)
goal = @spawn @time timeshift.get_field_shift(beh, spikes, c(-2:0.05:2); 
                                        multithread=true, 
                                        postfunc=info.information, newkws...);
goal_fineNarrow = @spawn @time timeshift.get_field_shift(beh, spikes, c(-1:0.01:1); 
                                        multithread=true, 
                                        postfunc=info.information, newkws...);
goal_broadTimes = @spawn @time timeshift.get_field_shift(beh, spikes, c(-4:0.05:4); 
                                        multithread=true, 
                                        postfunc=info.information, newkws...);
# C (( ---- SPECGOAL ---- ))
props = ["currentAngle", "currentPathLength", "stopWell"]
newkws = (; kws..., filters=merge(kws.filters, filters),
          resolution=[40, 40, 5], gaussian=2.3*0.5, props=props)
specgoal = @spawn @time timeshift.get_field_shift(beh, spikes, c(-1:0.05:1); 
                                        multithread=true, 
                                        postfunc=info.information, newkws...);
# C (( ---- SPECPLACE ---- ))
props = ["x", "y", "stopWell"]
newkws = (; kws..., filters=merge(kws.filters, filters),
          resolution=[40, 40, 5], gaussian=2.3*0.5, props=props)
specplace = @spawn @time timeshift.get_field_shift(beh, spikes, c(-1:0.05:1); 
                                        multithread=true, 
                                        postfunc=info.information, newkws...);
# D (( ---- FULL ---- ))
props = ["x", "y", "currentAngle", "currentPathLength", "stopWell"]
newkws = (; kws..., filters=merge(kws.filters, filters),
          resolution=[40, 40, 40, 40, 5], gaussian=2.3*0.5, props=props)
full = @spawn @time timeshift.get_field_shift(beh, spikes, c(-1:0.05:1); 
                                        multithread=true, 
                                        postfunc=info.information, newkws...);


# (( PLACE ))
if isdefined(Main, :place_broadTimes)
    place_broadTimes = Dict(fetch(place_broadTimes)...)
end
if isdefined(Main, :place_fineNarrow)
    place_fineNarrow = Dict(fetch(place_fineNarrow)...)
end
if isdefined(Main, :place)
    place = Dict(fetch(place)...)
end
if isdefined(Main, :place_correctOnly)
    place_correctOnly = Dict(fetch(place_correctOnly)...)
end

# (( GOAL ))
if isdefined(Main, :goal)
    goal = Dict(fetch(goal)...)
end
if isdefined(Main, :goal_fineNarrow)
    goal_fineNarrow = Dict(fetch(goal_fineNarrow)...)
end
if isdefined(Main, :goal_broadTimes)
    goal_broadTimes = Dict(fetch(goal_broadTimes)...)
end
# (( SPECGOAL ))
if isdefined(Main, :specgoal)
    specgoal = Dict(fetch(specgoal)...)
end
# (( SPECPLACE ))
if isdefined(Main, :specplace)
    specplace = Dict(fetch(specplace)...)
end
# (( FULL ))
if isdefined(Main, :full)
    full = Dict(fetch(full)...)
end
utils.pushover("Finished spawned shift processes")


plot_shifts(place_broadTimes,    desc="Traj length sample: ")
plot_shifts(place_fineNarrow,    desc="PLACE_HIGHRES")
plot_shifts(place,               desc="PLACE_LOWRES")
plot_shifts(goal_fineNarrow,     desc="GOAL_HIGHRES")
plot_shifts(goal_broadTimes,     desc="GOAL_TRAJ_LEN")
plot_shifts(specgoal,            desc="SPECGOAL")
plot_shifts(full,                desc="FULL")
plot_shifts(specplace,           desc="SPECPLACE")
plot_shifts(place_correctOnly,   desc="PLACE_CORRECTONLY")
plot_shifts(place_incorrectOnly, desc="PLACE_ERRORONLY")
plot_shifts(place_nontaskOnly,   desc="PLACE_NONTASK")
plot_shifts(place_correctCue,    desc="PLACE_CUEcorrect")
plot_shifts(place_correctMem,    desc="PLACE_MEMcorrect")

# Shuffle test
place = @time timeshift.get_field_shift_shufflesfield_shift(beh, spikes, -2:0.1:2; 
                                        multithread=true, 
                                        postfunc=info.information, newkws...);
