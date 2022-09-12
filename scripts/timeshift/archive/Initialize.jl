@time include(scriptsdir("fields", "Initialize.jl"))
info = Field.info
include(scriptsdir("timeshift", "TimeShift_setsOfInterest.jl"))
import Base.Threads: @spawn
using ThreadSafeDicts
using Blink
using Interact
using NaNStatistics
import .Timeshift
using .Timeshift.dataframe: info_to_dataframe
using Combinatorics: powerset
sp = copy(spikes)
convertToMinutes = true
_, spikes = Load.register(beh, spikes; transfer=["velVec"], on="time")

function get_key(;shifts, kws...)
    (;kws..., first=first(shifts), last=last(shifts), step=Float64(shifts.step)) 
end


# ----------
# Parameters
# ----------
PROPS = ["x", "y", "currentHeadEgoAngle", "currentPathLength", "stopWell"]
IDEALSIZE = Dict(key => (key=="stopWell" ? 5 : 40) for key in PROPS)
#-------------------------------------------------------
"""
Translate into shorcut names
"""
𝕄(items)  = [replace(item, recon_process.var_shortcut_names...)
             for item in items]
"""
UnTranslate from shorcut names
"""
𝕄̅(items)  = [replace(item, Dict(kv[2]=>kv[1] for kv in shortcut_names)...)
             for item in items]
sz(items) = [IDEALSIZE[item] for item in items]
#-------------------------------------------------------
splitby = ["unit", "area"]
prop_set = marginals_superhighprior_shuffle
filts = Filt.get_filters()
shifts = -4:0.2:4
shifts = convertToMinutes ? shifts./60 : shifts
