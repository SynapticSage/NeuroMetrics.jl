quickactivate("/home/ryoung/Projects/goal-code/")
@time include(scriptsdir("fields", "Include.jl"))
@time include(scriptsdir("fields", "Initialize.jl"))
includet(srcdir("shuffle.jl"))
includet(srcdir("field/info.jl"))
includet(srcdir("field/timeshift.jl"))
include(scriptsdir("fields", "TimeShift_setsOfInterest.jl"))
_, spikes = raw.register(beh, spikes; transfer=["velVec"], on="time")
import Base.Threads: @spawn
using ThreadSafeDicts
using .timeshift
using Combinatorics: powerset
sp = copy(spikes)
convertToMinutes = false

# ----------
# ----------
# ----------
# Parameters
# ----------
PROPS = ["x", "y", "currentHeadEgoAngle", "currentPathLength", "stopWell"]
IDEALSIZE = Dict(key => (key=="stopWell" ? 5 : 40) for key in PROPS)


#-------------------------------------------------------
"""
Translate into shorcut names
"""
𝕄(items)  = [replace(item, var_shortcut_names...) for item in items]
"""
UnTranslate from shorcut names
"""
𝕄̅(items)  = [replace(item, Dict(kv[2]=>kv[1] for kv in shortcut_names)...)
             for item in items]
sz(items) = [IDEALSIZE[item] for item in items]
#-------------------------------------------------------

splitby=["unit", "area"]
filters = merge(filt.notnan("currentAngle"), 
                filt.minmax("currentPathLength", 2, 150),
                filt.speed_lib, filt.cellcount)
gaussian=2.3*0.5

# COMPUTE INFORMATION @ DELAYS
I = Dict()
S = Dict()
for props ∈ marginals_highprior
    marginal = 𝕄(props)
    newkws = (; kws..., filters, splitby, gaussian, props,
              resolution=sz(props), multi=:single, postfunc=info.information)
    I[marginal] = @spawn @time timeshift.get_field_shift(beh, spikes, -2:0.1:2; 
                                                         newkws...)
    S[marginal] = @spawn @time timeshift.get_field_shift_shufflesfield_shift(beh, spikes, -2:0.1:2; 
                                                         newkws...)
end
for (key, value) in I
    I[key] = fetch(I[key])
end

# GET FIELDS AT BEST-(𝛕)
F = Dict()
for (key, value) in I
    newkws = (;kws..., filters, splitby, gaussian, props, resolution=sz(props))
    F      = timeshift.fetch_best_fields(I[key], beh, spikes; newkws...)
end


