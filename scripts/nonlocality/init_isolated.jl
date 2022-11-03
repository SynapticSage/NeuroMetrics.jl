(animal, day) = first((("RY22",21), ("super", 0)))

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

using DataStructures: OrderedDict
using DimensionalData
import DimensionalData: Between
using ProgressMeter
using DataFrames, DataFramesMeta
using Statistics, NaNStatistics, StatsBase, StatsPlots, HypothesisTests, GLM
using Plots
using LazySets


clab = OrderedDict(-1 => "nontask", 0 => "cue", 1=> "mem", missing=>"sleep")
Munge.nonlocal.setclab(clab)
isonames =  OrderedDict(false => :adjacent, true=>:isolated)
filt_desc = OrderedDict(:all => "> 2cm/s")
save_kws = (;pfc_rate_analy=true)
filt = Filt.get_filters()
datacut = :all

Plot.setappend((;animal,day))
Plot.setparentfolder("nonlocality")

# ===================
# ACQUIRE DATA
# ===================
# Acquire data
@time spikes, beh, cells = Load.load(animal, day, data_source=["spikes","behavior", "cells"])
GC.gc()
beh, spikes = Utils.filtreg.filterAndRegister(beh, spikes, on="time", transfer=["x","y","cuemem"], filters=filt[datacut], filter_skipmissingcols=true)
allspikes = copy(spikes)
beh2 = Load.load_behavior(animal,day)
Munge.nonlocal.setunfilteredbeh(beh2)

# Acquire LFP and isolated spikes
lfp = Load.load_lfp(animal, day, tet=:default);
lfp.time = lfp.time .- Load.min_time_records[1]
lfp = Munge.lfp.annotate_cycles(lfp)
@assert length(unique(lfp.cycle)) > 1
#sp = @subset(spikes, :tetrode .== 6);

F = load_fields()
kz = collect(filter(k->k.animal == animal && k.day == day, keys(F)))
@time f = F[bestpartialmatch(kz, 
                             (;datacut, widths=5,coactivity=nothing), 
                             nothing_means_removekey=true)];
f = f isa ShiftedFields ? f : ShiftedFields(deepcopy(f))
unitshift = Timeshift.types.matrixform(f)

# ===================
# OUT OF FIELD SPIKES
# ===================
# Setup a  shift-getting convenience method, the shifts, and a few metrics
shifts = collect(unitshift.dims[2])
push_dims!(unitshift)
push_celltable!( unitshift, cells, :unit, :area)
push_metric!(unitshift, Field.metrics.bitsperspike)
push_metric!(unitshift, Field.metrics.totalcount)
push_shiftmetric!(unitshift, best_tau!; metric=:bitsperspike)
annotate_nonlocal_spikes!(spikes, cells, unitshift, 0)

# ===================
# ISOLATED SPIKING
# ===================
#Munge.spiking.isolated(sp, lfp, include_samples=true)
Munge.spiking.isolated(spikes, lfp, include_samples=false)


# Which cells pass our criteria?
#region = :CA1
#metricfilter = metricfilters[region]
@assert length(unique(lfp.cycle)) > 1
cycles = Munge.lfp.get_cycle_table(lfp)
Munge.lfp.annotate_cycles!(spikes, cycles)

"""
Make sure sleep not included here
"""
tsk = Load.load_task(animal, day)
tsk = DataFrame(unique(eachrow(tsk[:,[:start,:end,:task]])))
tsk = subset(tsk, :task=>t->t .!= "sleep")
ismissing.(spikes.velVec)


using JLD2
filename = datadir("isolated","iso_animal=$(animal)_day=$(day)")
!isdir(dirname(filename)) ? mkpath(dirname(filename)) : nothing
@save "$filename"