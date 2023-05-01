begin
    using Markdown, InteractiveUtils, DrWatson, Revise
    quickactivate(expanduser("~/Projects/goal-code"))
    using DataFrames, DataFramesMeta, KernelDensity, Distributions, Plots,
          StatsPlots, Measures, Distributions, ProgressMeter, ProgressLogging,
          ThreadSafeDicts, NaNStatistics, Infiltrator, TimerOutputs, Serialization,
          DIutils
    using DataStructures: OrderedDict
    using Combinatorics: powerset
    import Base.Threads: @spawn
    using JLD2
    using GoalFetchAnalysis , 
           GoalFetchAnalysis.Timeshift, GoalFetchAnalysis.Timeshift.types, 
           GoalFetchAnalysis.Timeshift.checkpoint
    import DI: Filt
    using .Timeshift.dataframe: info_to_dataframe
    using .Field.recon_process: get_shortcutnames, inv_shortcutnames
end
opts = argparse()
keyfilterstr= k->!occursin("half",k) && !occursin("odds",k) # skip crossval
filts    = Filt.get_filters_precache(;keyfilterstr)
# datacuts = [:all, :cue, :memory, :task, :nontask, :cue_correct, :cue_error, :mem_correct,:mem_error, :correct,:error]
datacuts = collect(keys(filts))
sets = include(scriptsdir("timeshift","sets_of_interest.jl"))
prop_set = sets.marginals_superhighprior_shuffle
PROPS = ["x", "y", "currentHeadEgoAngle", "currentPathLength", "stopWell"]
IDEALSIZE = Dict(key => (key=="stopWell" ? 5 : 40) for key in PROPS)

WIDTHS = opts["width"]
thresh = opts["thresh"]
shifts= opts["init_shift"]:opts["period_shift"]:opts["final_shift"]
# datasets = (("RY16",36,nothing), ("RY22", 21, nothing)) #, ("RY16",36, :adj),("RY16",36, :iso),)
datasets = (("super_clean", 0, nothing),("super_clean",0,:iso), ("super_clean",0,:adj))
overwrite = false
isotet = :default
@assert datasets isa Tuple
function get_key(;shifts, kws...)
    (;kws..., grid=:adaptive,
    first=first(shifts), last=last(shifts), 
               step=Float64(shifts.step)) 
end
function keymessage(I::AbstractDict, key)
   @info key
   docontinue=false
   if DIutils.namedtup.orderlessmatch(key, keys(I))
       if I[key] isa Task && !(istaskfailed(I[key]))
           #@info "task key=$key already exists"
           printstyled("SKIPPING...\n", blink=true)
           docontinue=true
       elseif I[key] isa Task && istaskfailed(I[key])
           "key=$key already exists, but failed...redo!"
       else
           #@info "key=$key already exists"
           printstyled("SKIPPING...\n", blink=true)
           docontinue=true
       end
   end
   if key ∉ keys(I)
       #@info "key=$key ∉ keys, ...creating..."
   end
   docontinue
end
I = OrderedDict{NamedTuple, Any}()
F = OrderedDict{NamedTuple, Any}()
(animal, day, frac) = first(datasets)
global spikes = lfp = beh = cycles = nothing
(animal, day, frac) = first(datasets[1:end])
# ISSUE: MAKE SURE SUPER ANIMAL ON SOLID STATE
DIutils.pushover("Ready for cachemain.jl")

@showprogress "animal" for (animal, day, frac) in datasets[1:end] #DI.animal_set
    @info "loop" animal day frac
    @time spikes, beh, ripples, cells = DI.load(animal, day);
    if :animal ∉ propertynames(spikes)
        DIutils.filtreg.register(beh, spikes; transfer=["velVec"], on="time")
    else
        sort!(beh, [:animal,:time])
        sort!(spikes,[:animal,:time])
        DIutils.filtreg.register(groupby(beh, :animal),
                    groupby(spikes, :animal), transfer=["velVec"], on="time")
    end
    spikescopy = copy(spikes)
    if frac == :iso || frac == :adj
        lfp = DI.load_lfp(animal, day, tet=isotet, subtract_earlytime=true)
        @info "annotating cycles"
        @time Munge.lfp.annotate_cycles(lfp::DataFrame; 
                                  phase_col="phase", 
                                  method="peak-to-peak")
        # histogram(lfp.time); histogram!(spikes.time)
        @info "isolated cycles"
        spikes = Munge.spiking.isolated(spikes, lfp, refreshcyc=true)
        JLD2.@save "~/tmp_checkpoint_iso.jld2" {compress=true} spikes, lfp
        using Plot.lfplot
        @info "registering"
        _, spikes = DIutils.filtreg.register(lfp, spikes; transfer=["phase"], on="time")
        #phasepl
        @info "plotting"
        histogram(spikes.time, group=spikes.isolated, normalize=:pdf)
        histogram(spikes.phase, group=spikes.isolated, normalize=:pdf,
                  alpha=0.33)
        histogram(spikes.phase, group=spikes.isolated, normalize=:probability,
                  alpha=0.33)
        dropmissing!(spikes, :isolated)
        spikes = frac == :iso ? (println("selecting iso"); spikes[spikes.isolated, :]) : 
                                (println("selecting adj"); spikes[(!).(spikes.isolated), :])
        @info "unique iso" unique(spikes.isolated)
    end
    # CHECKPOINT
    shuffle_type = :dotson
    if shuffle_type == :dotson
        #nbins = 50
        Munge.behavior.annotate_relative_xtime!(beh)
        #beh.trajreltime_bin = floor.(beh.trajreltime * (nbins-1))
        if :animal ∉ propertynames(spikes)
            _, spikes = DIutils.filtreg.register(beh, spikes;
                                     transfer=["trajreltime","epoch"],
                                     on="time")
        elseif :animal ∈ propertynames(spikes)
            DIutils.filtreg.register(groupby(beh, :animal),
                                     groupby(spikes, :animal);
                                     transfer=["trajreltime","epoch"],
                                     on="time")
        end
        trajperiod = Table.get_periods(beh, :period)
        sort!(spikes,:time)
        spike_trajs = DIutils.searchsortedprevious.([trajperiod.start], spikes.time)
        spikes.trajstart, spikes.trajdel = eachcol(trajperiod[spike_trajs,[:start, :δ]])
        mean.(map(x-> x.>0 .* x .< spikes.trajdel, [spikes.time .- spikes.trajstart]))
    end
    spikes, beh = copy(spikes), copy(beh); GC.gc() 
    datacuts = collect(keys(filts))
    # LOAD PREXISTING?
    #if isfile(Timeshift.mainspath())
    #    I = Timeshift.load_mains()
    #    F = OrderedDict{NamedTuple, Any}()
    #else
    #end
    #keys(I)
    statfile = datadir("exp_pro", "timeshift", "cachelog.txt")
    key=(;)
    result_dict = OrderedDict{Real, Any}()
    spikes = dropmissing(spikes, :trajreltime)
    iter = Iterators.product(WIDTHS, thresh, datacuts, prop_set)
    sort!(beh, [:time]) # required for multi-animal datasets
    sort!(spikes, [:time])
    @showprogress "cuts and params" for (widths, thresh, datacut, props) ∈ iter
        try
            marginal = get_shortcutnames(props)
            key      = get_key(;marginal, datacut, shifts, widths, thresh,
                                animal, day, frac)
            filt = filts[datacut]
            @info filt filts[datacut]
            if keymessage(I, key); continue; end
            @time tmp = Timeshift.shifted_fields(beh, spikes,
                        shifts, props; fieldpreset=:yartsev,
                                       shiftbeh=false,
                                       widths, filters=filt, thresh)
            tmp = Timeshift.DictOfShiftOfUnit{ keytype(tmp)}(tmp)
            F[key] = ShiftedFields(tmp)
            I[key] = F[key].metrics
        catch exception
            # @infiltrate
            @info "exception" exception
        end
    end
    DIutils.pushover("Finished cachemain.jl $animal $day $frac",title="Timeshift")
    # @infiltrate
    # SAVE the shifting information
    # archivestr = isempty(setdiff(unique([d[3] for d in datasets]), 
    # [:adj,:iso])) ? "iso" : ""
    archivestr = frac !== nothing ? [string(animal), string(day), string(frac), string(isotet)] :
                    [string(animal), string(day)]
    archivestr  = join(archivestr, "-")
    checkpoint.save_fields(F; overwrite, archive=archivestr);
end

# Make a dataframe of the namedtuple keys of F
keymat = hcat([(key|>collect) for key in keys(F)]...)
keymat = permutedims(keymat, (2,1))
keyframe = DataFrame(keymat,
    keys(F)|>first|>propertynames.|>string|>collect)

# Let's make a list of properties that change
changing = [name for name in names(keyframe) 
    if length(unique(keyframe[!,name])) > 1]
# And add some things I might change over multiple iterations
changing = [changing..., "animal", "step"]

# Now we want to extract a column of timeshifts for every single set
# of processed data, provid a column name based on information in `changing`,
# then we will want to 
using GoalFetchAnalysis.Munge.timeshift
(key, value) = F |> first;
for (key, value) in F
    savekey = OrderedDict(k=>v for (k,v) in Dict(pairs(key)) 
                      if k in Symbol.(changing))
    str   = DIutils.dict.tostring(savekey)
    value = ShiftedFields(value)
    value = matrixform(value)
    # BUG: `unit` values associated with each `ShiftedField` object is wrong...probably
    V = timeshift.usefulmetrics(value, cells)
end

# In this section, we optionall add information to the cells table
# about the time 
@time F = checkpoint.load_fields()
map(keys(F)|>collect) do k
    k.animal
end |> unique
map(keys(F)|>collect) do k
    k.frac
end |> unique
keys(F)|>collect

DI.superanimal.superanimal_clean_times_and_neurons(0; append="_clean")
printstyled("FInished superanimal_clean_times_and_neurons", blink=true)

# Convert F to dataframe
@time Fd = Table.to_dataframe(F,       key_name=["shift"], explode=false);
Fd
