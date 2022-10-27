using Markdown
using InteractiveUtils

using DrWatson, Revise
quickactivate(expanduser("~/Projects/goal-code"))
using PlutoUII
PlutoUI.TableOfContents(title="Caching Mains and Shuffles")

using DataFrames, DataFramesMeta
using DataStructures: OrderedDict
using KernelDensity, Distributions
using Plots, StatsPlots, Measures, Distributions
using ProgressMeter, ProgressLogging
using Combinatorics: powerset
import Base.Threads: @spawn
using ThreadSafeDicts, NaNStatistics
PlutoUI.TableOfContents(title="Caching Mains and Shuffles")

using DataFrames, DataFramesMeta
using DataStructures: OrderedDict
using KernelDensity, Distributions
using Plots, StatsPlots, Measures, Distributions
using ProgressMeter, ProgressLogging
using Combinatorics: powerset
import Base.Threads: @spawn
using ThreadSafeDicts, NaNStatistics
using Infiltrator
using TimerOutputs

using GoalFetchAnalysis 
import Timeshift
using Timeshift.dataframe: info_to_dataframe
using Field.recon_process: get_shortcutnames, inv_shortcutnames
import Load
import Filt
import Munge
filts = Filt.get_filters_precache()

# source:  https://github.com/fonsp/Pluto.jl/issues/115  
function ingredients(path::String)
    # this is from the Julia source code (evalfile in base/loading.jl)
    # but with the modification that it returns the module instead of the last object
    name = Symbol(basename(path))
    m = Module(name)
    Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
    m
end
sets = ingredients(scriptsdir("timeshift", "TimeShift_setsOfInterest.jl"))
prop_set = sets.marginals_superhighprior_shuffle


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
 PROPS = ["x", "y", "currentHeadEgoAngle", "currentPathLength", "stopWell"]
 IDEALSIZE = Dict(key => (key=="stopWell" ? 5 : 40) for key in PROPS)
 shifts=-2:0.05:2
 widths = 5.0f0
 thresh = 1.5f0

shuffle_type = :dotson
@time spikes, beh, ripples, cells = Load.load("RY16", 36);
_, spikes = Load.register(beh, spikes; transfer=["velVec"], on="time")
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

#if isfile(Timeshift.shufflespath())
#    S = Timeshift.load_shuffles()
#else
#    S = OrderedDict{NamedTuple, Any}()
#end
#keys(S)

result_dict = nothing
spikes         = dropmissing(spikes, :trajreltime)

for datacut in [:memory, :cue]

    props = first( prop_set)
    marginal = get_shortcutnames(props)
    key = get_key(;marginal, datacut, shifts, widths, thresh)
    filt = filts[datacut]

    result_dict = OrderedDict{NamedTuple, 
                              Timeshift.DictOfShiftOfUnit{Float64}}()
    @showprogress for shuffle in 1:100

        Timeshift.shuffle.shifted_field_shuffles(beh, spikes,
                    shifts, props; 
                    fieldpreset=:yartsev,
                    shufflepreset=shuffle_type, 
                    nShuffle=1,
                    savename="datacut=$datacut",
                    saveAfter=1,
                    widths, 
                    thresh,
                    shiftbeh=false,
                    filters=filt,
                    result_dict)
    end
end
