quickactivate("/home/ryoung/Projects/goal-code/")
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
convertToMinutes = true

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
𝕄(items)  = [replace(item, recon_process.var_shortcut_names...) 
             for item in items]
"""
UnTranslate from shorcut names
"""
𝕄̅(items)  = [replace(item, Dict(kv[2]=>kv[1] for kv in shortcut_names)...)
             for item in items]
sz(items) = [IDEALSIZE[item] for item in items]
#-------------------------------------------------------

splitby=["unit", "area"]
#gaussian=2.3*0.5
filts = filt.get_filters()
shifts = -4:0.2:4
shifts = convertToMinutes ? shifts./60 : shifts
function get_key(;shifts, kws...)
    (;kws..., first=first(shifts), last=last(shifts), step=Float64(shifts.step)) 
end

# List of the ways that one might want to vary this analysis
# 1. filters
#   all-times
#       task
#           correct
#           error
#           cue -> correct/error
#           mem -> correct/error
#       nontask
# 2. marginals
#   ["x", "y"]
#   ["currentHeadEgoAngle", "currentPathLength"]
#   ["currentHeadEgoAngle", "currentPathLength", "stopWell"]
#   ["x", "y", "currentHeadEgoAngle", "currentPathLength", "stopWell"]
# 3. Resolution
#   -8 : 8 seconds (trajectory length)
#   (block length)

I = OrderedDict()
# COMPUTE INFORMATION @ DELAYS
@showprogress 0.1 "Datacut iteration" for datacut ∈ keys(filts)
    for props ∈ marginals_highprior
        marginal = 𝕄(props)
        key = get_key(;marginal, datacut, shuf=:cDt_t, shifts)
        if key in keys(I)
            if I[key] isa Task && !(istaskfailed(I[key]))
                @info "key=$key already exists, skipping..."
                continue
            else
                @info "key=$key already exists, but failed...redo!"
            end
        end
        newkws = (; kws..., filters=filts[datacut], splitby, props, dokde=false,
                  resolution=sz(props), multi=:single, postfunc=info.information)
        I[key] = @spawn @time timeshift.get_field_shift(beh, spikes, shifts; 
                                                             newkws...)
    end
    try
        for (key, value) in I
            @info I
            I[key] = fetch(I[key])
        end
        @info I
        utils.pushover("Done fetchging jobs for datacut=$datacut")
        timeshift.save_mains(I)
    catch
        @warn "potential task failure for props=$props, datacut=$datacut"
    end
end


# COMPUTE SHUFFLE INFORMATION @ DELAYS
S = OrderedDict()

@showprogress 0.1 "Datacut shuffle iteration" for datacut ∈ keys(filts)
    for props ∈ marginals_highprior
        marginal = 𝕄(props)
        key = get_key(;marginal, datacut, shifts, shuffle_type="σ(traj)")
        if key in keys(S)
            if S[key] isa Task && !(istaskfailed(S[key]))
                @info "key=$key already exists, skipping..."
                continue
            else
                @info "key=$key already exists, but failed...redo!"
            end
        end
        newkws = (; kws..., filters=filts[datacut], splitby, props,
                  resolution=sz(props), multi=:single,
                  exfiltrateAfter=20,
                  postfunc=info.information)
        S[key] = @time timeshift.get_field_shift_shuffles(beh, spikes, shifts; 
                                                         newkws...)
    end

    try
        for (key, value) in S
            @info S
            S[key] = fetch(S[key])
        end
        @info S
        utils.pushover("Done fetchging jobs for datacut=$datacut")
    catch
        @warn "potential task failure for props=$props, datacut=$datacut"
    end
    timeshift.save_shuffles(S)

end

# GET FIELDS AT BEST-(𝛕)
F = Dict()
for (key, value) in I
    newkws = (;kws..., filters=filts, splitby, gaussian, props, resolution=sz(props))
    F      = timeshift.fetch_best_fields(I[key], beh, spikes; newkws...)
end


