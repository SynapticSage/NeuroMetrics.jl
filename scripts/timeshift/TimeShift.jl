quickactivate("/home/ryoung/Projects/goal-code/")
@time include(scriptsdir("fields", "Initialize.jl"))
includet(srcdir("shuffle.jl"))
info = field.info
include(scriptsdir("timeshift", "TimeShift_setsOfInterest.jl"))
_, spikes = raw.register(beh, spikes; transfer=["velVec"], on="time")
import Base.Threads: @spawn
using ThreadSafeDicts
using .timeshift
using Combinatorics: powerset
sp = copy(spikes)
convertToMinutes = true

# ----------
# Parameters
# ----------
PROPS = ["x", "y", "currentHeadEgoAngle", "currentPathLength", "stopWell"]
IDEALSIZE = Dict(key => (key=="stopWell" ? 5 : 40) for key in PROPS)
#-------------------------------------------------------
"""
Translate into shorcut names
"""
𝕄(items)  = [replace(item, recon_process.var_shortcut_names...) for item in items]
"""
UnTranslate from shorcut names
"""
𝕄̅(items)  = [replace(item, Dict(kv[2]=>kv[1] for kv in shortcut_names)...) for item in items]
sz(items) = [IDEALSIZE[item] for item in items]
#-------------------------------------------------------
splitby=["unit", "area"]
prop_set = marginals_superhighprior_shuffle
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


# ========
# MAINS
# ========
# Main, non-shuffles
if isfile(timeshift.mainspath())
    I = timeshift.load_mains()
else
    I = OrderedDict()
end
# COMPUTE INFORMATION @ DELAYS
@showprogress 0.1 "Datacut iteration" for datacut ∈ keys(filts)
    running = false
    @info "Datacut = $datacut"
    for props ∈ 
        @info "Props = $props"
        marginal = 𝕄(props)
        key = get_key(;marginal, datacut, shifts)
        if key ∈ keys(I)
            if I[key] isa Task && !(istaskfailed(I[key]))
                @info "task key=$key already exists"
                printstyled("SKIPPING...\n", blink=true)
                continue
            elseif I[key] isa Task && istaskfailed(I[key])
                "key=$key already exists, but failed...redo!"
            else
                @info "key=$key already exists"
                printstyled("SKIPPING...\n", blink=true)
                continue
            end
        end
        if key ∉ keys(I)
            @info "key=$key ∉ keys, ...creating..."
        end
        newkws = (; kws..., filters=filts[datacut], splitby, props, dokde=false,
                  resolution=sz(props), multi=:single, postfunc=info.information)
        try
            I[key] = @spawn @time timeshift.get_field_shift(beh, spikes, shifts; newkws...)
            I[key] = fetch(I[key])
        catch
            @warn "key=$key does not run, possibly an edge case where your filters are too stringent for the behavioral property measured"
        end
        running = true
    end
    if running
        try
            for (key, value) in I
                @info I
                I[key] = fetch(I[key])
            end
            @info I
            utils.pushover("Done fetchging jobs for datacut=$datacut")
            timeshift.save_mains(I)
            @infiltrate
        catch
            @warn "potential task failure for props=$props, datacut=$datacut"
        end
    end
end

# ========
# SHUFFLE
# ========
shuffle_type = :dotson
if shuffle_type == :dotson
    nbins = 50
    raw.behavior.annotate_relative_xtime!(beh)
    beh.trajreltime_bin = floor.(beh.trajreltime * (nbins-1))
    _, spikes = raw.register(beh, spikes; 
                             transfer=["trajreltime","trajreltime_bin"], 
                             on="time")
end
# COMPUTE SHUFFLE INFORMATION @ DELAYS
if isfile(timeshift.shufflespath())
    S = timeshift.load_shuffles()
else
    S = OrderedDict()
end
@showprogress 0.1 "Datacut shuffle iteration" for datacut ∈ keys(filts)
    for props ∈ prop_set
        marginal = 𝕄(props)
        key = get_key(;marginal, datacut, shifts, shuffle_type)
        if key in keys(S)
            if (S[key] isa Task && !(istaskfailed(S[key]))) ||
                !(S[key] isa Task)
                @info "key=$key already exists, skipping..."
                continue
            else
                @info "key=$key already exists, but failed...redo!"
            end
        else
                @info "key=$key"
        end
        shufshift_settings = (; kws..., 
                    filters=filts[datacut], 
                    preset=shuffle_type,
                    splitby, props,
                    resolution=sz(props), 
                    compute=:single,
                    exfiltrateAfter=100,
                    postfunc=info.information)
        S[key] = @time timeshift.get_field_shift_shuffles(beh, spikes, shifts; shufshift_settings...)
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

#

# ================
# Get fields at shifts
# ================
# Main, non-shuffles
if isfile(timeshift.fieldspath())
    F = timeshift.load_fields()
else
    F = OrderedDict()
end

# COMPUTE INFORMATION @ DELAYS
@showprogress 0.1 "Datacut iteration" for datacut ∈ keys(filts)
    running = false
    @info "Datacut = $datacut"
    for props ∈ prop_set
        @info "Props = $props"
        marginal = 𝕄(props)
        key = get_key(;marginal, datacut, shifts)
        if key ∈ keys(F)
            if F[key] isa Task && !(istaskfailed(I[key]))
                @info "task key=$key already exists"
                printstyled("SKIPPING...\n", blink=true)
                continue
            elseif F[key] isa Task && istaskfailed(I[key])
                "key=$key already exists, but failed...redo!"
            else
                @info "key=$key already exists"
                printstyled("SKIPPING...\n", blink=true)
                continue
            end
        end
        if key ∉ keys(F)
            @info "key=$key ∉ keys, ...creating..."
        end
        newkws = (; kws..., filters=filts[datacut], splitby, props, dokde=false,
                  resolution=sz(props), multi=:single)
        try
            F[key] = @time timeshift.get_field_shift(beh, spikes, shifts; newkws...)
            F[key] = fetch(I[key])
        catch
            @warn "key=$key does not run, possibly an edge case where your filters are too stringent for the behavioral property measured"
        end
        running = true
    end
    if running
        try
            for (key, value) in F
                @info F
                F[key] = fetch(I[key])
            end
            @info F
            utils.pushover("Done fetchging jobs for datacut=$datacut")
            timeshift.save_fields(F)
        catch
            @warn "potential task failure for props=$props, datacut=$datacut"
        end
    end
end

# ================
# Examine best TAU
# ================

# Place best tau into cell.csv storage
key = (;first(keys(I))..., datacut=:all)
cellTaus = sort(timeshift.imax(info_to_dataframe(I[key], shift_scale=:minuntes)),:unit)[:,[:unit,:area,:bestTau]]
keyname = keytostring(key, ", "=>"_"; remove_numerics=true)
rename!(cellTaus, "bestTau" => "bestTau_$keyname")
raw.save_cell_taginfo(cellTaus, animal, day, keyname)


# GET FIELDS AT BEST-(𝛕)
F = Dict()
for (key, value) in I
    newkws = (;kws..., filters=filts, splitby, gaussian, props, resolution=sz(props))
    F      = timeshift.fetch_best_fields(I[key], beh, spikes; newkws...)
end


# ==============
# Cross-validate
# ==============

