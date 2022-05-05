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


jobs = Dict()
n = 8

# A (( ---- PLACE ---- ))
props = ["x", "y"]
splitby=["unit", "area"]
kws=(;splitby, filters=merge(filt.speed_lib, filt.cellcount), multi=:single, postfunc=info.information)

newkws = (; kws..., resolution=40, gaussian=2.3*0.5, props=props, filters=merge(kws.filters))
jobs[:place] = @spawn @time timeshift.get_field_shift(beh, spikes, c(-n:0.1:n); newkws...);

newkws = (; kws..., resolution=40, gaussian=4.3*0.5, props=props, filters=merge(kws.filters, filt.correct))
jobs[:place_correctOnly] = @time timeshift.get_field_shift(beh, spikes, c(-n:0.1:n); newkws...);

newkws = (; kws..., resolution=40, gaussian=4.3*0.5, props=props, filters=merge(kws.filters, filt.incorrect))
jobs[:place_incorrectOnly] = @spawn @time timeshift.get_field_shift(beh, spikes, c(-n:0.1:n); newkws...);

newkws = (; kws..., resolution=40, gaussian=4.3*0.5, props=props, filters=merge(kws.filters, filt.nontask))
jobs[:place_nontaskOnly] = @spawn @time timeshift.get_field_shift(beh, spikes, c(-n:0.1:n); newkws...);

newkws = (; kws..., resolution=40, gaussian=4.3*0.5, props=props, filters=merge(kws.filters, filt.cue, filt.correct))
jobs[:place_correctCue] = @spawn @time timeshift.get_field_shift(beh, spikes, c(-n:0.1:n); newkws...);

newkws = (; kws..., resolution=40, gaussian=4.3*0.5, props=props, filters=merge(kws.filters, filt.mem, filt.correct))
jobs[:place_correctMem] = @spawn @time timeshift.get_field_shift(beh, spikes, c(-n:0.1:n); newkws...);


# B (( ---- GOAL ---- ))
filters = merge(kws.filters, filt.correct,
                filt.notnan("currentAngle"), 
                filt.minmax("currentPathLength", 4, 150))
props = ["currentAngle", "currentPathLength"]

newkws = (; kws..., filters=merge(kws.filters, filters), resolution=40, gaussian=4.3*0.5, props=props)
jobs[:goal] = @spawn @time timeshift.get_field_shift(beh, spikes, c(-n:0.1:n); newkws...);

# C (( ---- SPECGOAL ---- ))
props = ["currentAngle", "currentPathLength", "stopWell"]
newkws = (; kws..., filters=merge(kws.filters, filters),
          resolution=[40, 40, 5], gaussian=n.3*0.5, props=props)
jobs[:specgoal] = @spawn @time timeshift.get_field_shift(beh, spikes, c(-n:0.1:n); newkws...);
# C (( ---- SPECPLACE ---- ))
props = ["x", "y", "stopWell"]
newkws = (; kws..., filters=merge(kws.filters, filters),
          resolution=[40, 40, 5], gaussian=n.3*0.5, props=props)
jobs[:specplace] = @spawn @time timeshift.get_field_shift(beh, spikes, c(-n:0.1:n); 
                                        newkws...);
# D (( ---- FULL ---- ))
props = ["x", "y", "currentAngle", "currentPathLength", "stopWell"]
newkws = (; kws..., filters=merge(kws.filters, filters),
          resolution=[40, 40, 40, 40, 5], gaussian=n.3*0.5, props=props)
jobs[:full] = @spawn @time timeshift.get_field_shift(beh, spikes, c(-n:0.1:n); newkws...);


# Retrieve
for (key,job) in jobs
    jobs[key] = Dict(fetch(job)...)
    #utils.pushover("Finished $key")
end


for (key,job) in jobs
    plot_shifts(job; desc=String(key), clim=:cell)
end

# Shuffle test
place = @time timeshift.get_field_shift_shufflesfield_shift(beh, spikes, -2:0.1:2; 
                                        multithread=true, 
                                        postfunc=info.information, newkws...);
