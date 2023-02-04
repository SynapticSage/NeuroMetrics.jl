
using Markdown, InteractiveUtils, DrWatson, Revise
quickactivate(expanduser("~/Projects/goal-code"))

using DataFrames, DataFramesMeta, KernelDensity, Distributions, Plots,
      StatsPlots, Measures, Distributions, ProgressMeter, ProgressLogging,
      ThreadSafeDicts, NaNStatistics, Infiltrator, TimerOutputs, Serialization

using DataStructures: OrderedDict
using Combinatorics: powerset
import Base.Threads: @spawn

using GoalFetchAnalysis , Timeshift, Timeshift.types
import Load, Filt, Munge
using Timeshift.dataframe: info_to_dataframe
using Field.recon_process: get_shortcutnames, inv_shortcutnames

opts = argparse()

keyfilterstr= k->!occursin("half",k) && !occursin("odds",k) # skip crossval
filts    = Filt.get_filters_precache(;keyfilterstr)
# datacuts = [:all, :cue, :memory, :task, :nontask, :cue_correct, :cue_error, :mem_correct,:mem_error, :correct,:error]
datacuts = collect(keys(filts))

sets = include(scriptsdir("timeshift","sets_of_interest.jl"))
prop_set = sets.marginals_superhighprior_shuffle

PROPS = ["x", "y", "currentHeadEgoAngle", "currentPathLength", "stopWell"]
IDEALSIZE = Dict(key => (key=="stopWell" ? 5 : 40) for key in PROPS)

WIDTHS = opt["width"]
thresh = opt["thresh"]
shifts= opt["init_shift"]:opt["period_shift"]:opt["final_shift"]
datasets = (("RY16",36,nothing), ("RY22", 21, nothing)) #, ("RY16",36, :adj),("RY16",36, :iso),)

function get_key(;shifts, kws...)
    (;kws..., grid=:adaptive,
    first=first(shifts), last=last(shifts), 
               step=Float64(shifts.step)) 
end
function keymessage(I::AbstractDict, key)
   @info key
   docontinue=false
   if Utils.namedtup.orderlessmatch(key, keys(I))
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
#

I = OrderedDict{NamedTuple, Any}()
F = OrderedDict{NamedTuple, Any}()
(animal, day, frac) = first(datasets)

@showprogress "animal" for (animal, day, frac) in datasets #Load.animal_set

    @info "loop" animal day frac

    @time spikes, beh, ripples, cells = Load.load(animal, day);
    _, spikes = Load.register(beh, spikes; transfer=["velVec"], on="time")


    if frac == :iso || frac == :adj
        lfp = Load.load_lfp(animal, day, tet=:ca1ref, subtract_earlytime=true)
        Munge.lfp.annotate_cycles(lfp::DataFrame; 
                                  phase_col="phase", 
                                  method="peak-to-peak")
        histogram(lfp.time); histogram!(spikes.time)
        spikes = Munge.spiking.isolated(spikes, lfp, refreshcyc=true)
        using Plot.lfplot
        _, spikes = Load.register(lfp, spikes; transfer=["phase"], on="time")
        #phasepl
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

    shuffle_type = :dotson
    if shuffle_type == :dotson
        #nbins = 50
        Munge.behavior.annotate_relative_xtime!(beh)
        #beh.trajreltime_bin = floor.(beh.trajreltime * (nbins-1))
        _, spikes = Load.register(beh, spikes;
                                 transfer=["trajreltime","epoch"],
                                 on="time")
        trajperiod = Table.get_periods(beh, :period)
        spike_trajs = Utils.searchsortedprevious.([trajperiod.start], spikes.time)
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
    @showprogress "widths" for widths ∈ WIDTHS
        @showprogress "Datacut iteration" for datacut ∈ datacuts
            #[:all, :cue, :memory, :task, :nontask]
            finished_batch = false
            @showprogress "Props" for props ∈ prop_set
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
                    @info "exception" exception
                end
            end
        end
        #Timeshift.save_mains(I)
        #Timeshift.save_mains(F)
    end
end

#savefile = datadir("timeshift","fixed_shifts_$shifts.serial")
#serialize(savefile, (;F,I,shifts))
#(F,I,shifts) = deserialize(savefile);

overwrite = true
archive = isempty(setdiff(unique([d[3] for d in datasets]), [:adj,:iso])) ? "iso" : ""
checkpoint.save_fields(F; overwrite, archive);

