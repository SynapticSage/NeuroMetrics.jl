
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

using JLD2
filename = datadir("isolated","iso_animal=$(animal)_day=$(day)")
load(filename)

for (k,v) in data
   println(k)
   Core.eval(Main, :($(Symbol(k)) = $v))
end

