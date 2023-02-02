@time include(scriptsdir("timeshift", "Initialize.jl"))
import Timeshift

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
if isfile(Timeshift.mainspath())
    I = Timeshift.load_mains()
else
    I = OrderedDict()
end

# COMPUTE INFORMATION @ DELAYS
@showprogress 0.1 "Datacut iteration" for datacut ∈ keys(filts)
    running = false
    @info "Datacut = $datacut"
    for props ∈ prop_set
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
        #try
            I[key] = @time Timeshift.get_field_shift(beh, spikes, shifts; newkws...)
            I[key] = fetch(I[key])
        #catch
            #@warn "key=$key does not run, possibly an edge case where your filters are too stringent for the behavioral property measured"
        #end
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
            Timeshift.save_mains(I)
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
    Munge.behavior.annotate_relative_xtime!(beh)
    beh.trajreltime_bin = floor.(beh.trajreltime * (nbins-1))
    _, spikes = Load.register(beh, spikes; 
                             transfer=["trajreltime","trajreltime_bin"], 
                             on="time")
end

# COMPUTE SHUFFLE INFORMATION @ DELAYS
if isfile(Timeshift.shufflespath())
    S = Timeshift.load_shuffles()
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
                println()
                @info "key=$key already exists, skipping..."
                println()
                continue
            else
                println()
                @info "key=$key already exists, but failed...redo!"
                println()
            end
        else
                println()
                @info "key=$key"
                println()
        end
        shufshift_settings = (; kws..., 
                    filters=filts[datacut], 
                    preset=shuffle_type,
                    splitby, props,
                    resolution=sz(props), 
                    compute=:single,
                    exfiltrateAfter=100,
                    postfunc=info.information)
        S[key] = @time Timeshift.get_field_shift_shuffles(beh, spikes, shifts;
                                                          shufshift_settings...)
    end

    try
        for (key, value) in S
            @info S
            S[key] = fetch(S[key])
        end
        @info S
        #utils.pushover("Done fetchging jobs for datacut=$datacut")
    catch
        @warn "potential task failure for props=$props, datacut=$datacut"
    end
    Timeshift.save_shuffles(S)
end

# ================
# FIELDS
# ================
# Main, non-shuffles
if isfile(Timeshift.fieldspath())
    F = Timeshift.load_fields()
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
            F[key] = @time Timeshift.get_field_shift(beh, spikes, shifts; newkws...)
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
            Timeshift.save_fields(F)
        catch
            @warn "potential task failure for props=$props, datacut=$datacut"
        end
    end
end

# Posthoc adding labels to namped-tuples
I, S, F = Timeshift.load_mains(),
          Timeshift.load_shuffles(),
          Timeshift.load_fields()
          
getfirstkey(x) = x[first(keys(x))]
firstkey(x)    = first(keys(x))
function add_namedtuple_fields(X::AbstractDict, detail::NamedTuple)
    for key in keys(X)
        x = pop!(X, key)
        key = (;key..., detail...)
        X[key] = x
    end
end

detail = (;grid=:fixed, resolution=40)
add_namedtuple_fields(I, detail)
add_namedtuple_fields(S, detail)
add_namedtuple_fields(F, detail)

Timeshift.save_mains(I, overwrite=true)
Timeshift.save_shuffles(S, overwrite=true)
Timeshift.save_fields(F, overwrite=true)

# Compare FIELDS AT BEST-(𝛕)
ui = @manipulate for i ∈ 1:size(imax,1)
    row = imax[i,:]
    fr = round(cells[row.unit,:meanrate];digits=2)
    zerobit = @subset(iall, :shift .==0 .&& :unit .== row.unit).info[1]
    bestbit = row.info_maximum
    T = Plots.plot(;plot_title="cell=$(row.unit), μ(fr)=$fr", 
                   grid=false, showaxis=false, yticks=[],
                   xticks=[], bottom_margin = -50Plots.px)
    layout = Plots.@layout [A{0.2h}; a b]
    best = @subset(f_select, :shift.==row.bestTauOrig .&& :unit .== row.unit)[1,:]
    zero = @subset(f_select, :shift.==0 .&& :unit .== row.unit)[1,:]
    clim_0   = ([nanquantile(vec(zero.value),x) for x in (0.20, 0.99)]...,)
    clim_max = ([nanquantile(vec(best.value),x) for x in (0.20, 0.99)]...,)
    bestheat = heatmap(best.value, title="best, bits=$bestbit", clim=clim_max)
    zeroheat = heatmap(zero.value, title="0, bits=$zerobit", clim=clim_0)
    plot(T,
        bestheat,
        zeroheat,
        layout=layout,
        size=(1400,900)
    )
end
w=Window()
body!(w,ui)

cell_plots=[]
for i = 1:size(imax,1)
    cell_plot=[]
    for shift in sort(unique(imax.bestTauOrig))
        hm = nothing
        try
            hm = @subset(f_select, :units .== imax[i,:unit] .&& :shift .== shift)[1,:].value
        catch
            @infiltrate
            @warn "Failed"
            hm = heatmap(NaN*ones(1,1))
        end
        push!(cell_plot, heatmap(hm, title="shift=$shift", xticks=[], yticks=[],
                                 colorbar=:none))
    end
    push!(cell_plots, cell_plot)
end

ui = @manipulate for i = 1:size(cells,1)
    p = try
        plot(cell_plots[i]..., size=(1000,1000))
    catch
        @warn "failed"
        plot()
    end
end
w=Window()
body!(w,ui)

# ==============
# Cross-validate
# ==============

# ================
# Examine best TAU
# ================

# Place best tau into cell.csv storage
key_order = (;marginal=key.marginal, 
             datacut=key.datacut,
             first=key.first,
             last=key.last,
             step=key.step)
cellTaus = sort(Timeshift.imax(info_to_dataframe(I[key], shift_scale=:minutes)),
                :unit)[:,[:unit,:area,:bestTau]]
keyname = Timeshift.keytostring(key, ", "=>"_"; remove_numerics=true)
rename!(cellTaus, "bestTau" => "bestTau_$keyname")
Load.save_cell_taginfo(cellTaus, animal, day, keyname)

utils.bestpartialmatch(F, key)
f_select  = table.to_dataframe(utils.applyvalues(, x->x.Cₕ);
                               key_name="shift", explode=false)
