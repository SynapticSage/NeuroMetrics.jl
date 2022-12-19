quickactivate(expanduser("~/Projects/goal-code/"));
using GoalFetchAnalysis
using Timeshift
using Timeshift.types
using Timeshift.shiftmetrics
using Field.metrics
using Plot
using Plot.receptivefield
using Utils.namedtup
using Munge.timeshift: getshift
using Munge.nonlocal
using Utils.statistic: pfunc
import Plot
using Munge.spiking
using Filt
using Infiltrator
using SoftGlobalScope

using DataStructures: OrderedDict
using DimensionalData
import DimensionalData: Between
using ProgressMeter
using DataFrames, DataFramesMeta
using Statistics, NaNStatistics, StatsBase, StatsPlots, HypothesisTests, GLM
using Plots
using LazySets
GC.gc()
datasets = (("RY16",36, nothing),("RY22",21, ""), ("RY22",21, "ref"))
for (animal,day,ref) in datasets[2:end]

    clab = OrderedDict(-1 => "nontask", 0 => "cue", 1=> "mem", missing=>"sleep")
    Munge.nonlocal.setclab(clab)
    isonames =  OrderedDict(false => :adjacent, true=>:isolated)
    filt_desc = OrderedDict(:all => "> 2cm/s")
    save_kws = (;pfc_rate_analy=true)
    filt = Filt.get_filters()
    datacut = :all

    Plot.setappend((;animal,day,ref))
    Plot.setparentfolder("nonlocality")

    # ===================
    # ACQUIRE DATA
    # ===================
    # Acquire data
    @time spikes, beh, cells = Load.load(animal, day, data_source=["spikes","behavior", "cells"])
    GC.gc()
    beh, spikes = Utils.filtreg.filterAndRegister(beh, spikes, on="time", transfer=["x","y","cuemem"], 
    filters=filt[datacut], filter_skipmissingcols=true)
    allspikes = copy(spikes)
    beh2 = Load.load_behavior(animal,day)
    Munge.nonlocal.setunfilteredbeh(beh2)

    # THis section requires enormous RAM, 100GB at least
    LFP = Load.load_lfp(animal, day; ref)
    ca1 = sort(unique(@subset(cells, :area .== "CA1").tetrode))
    LFP = groupby(LFP, :tetrode)
    LFP = [LFP[(;tetrode)] for tetrode in ca1]
    GC.gc()
    LFP = vcat(LFP...)
    GC.gc()

    Load.lfp.split_lfp_by_tet(animal, day; lfp=LFP, ref)
    LFP = nothing
    GC.gc()
end
