using GoalFetchAnalysis
using .Timeshift, .Plot, .Timeshift.types, .Timeshift.shiftmetrics, 
      .Field.metrics, .Plot.receptivefield, .DIutils.namedtup, 
      .Munge.isolated, .Munge.nonlocal, .Munge.spiking, .Plot.lfplot,
      .DIutils.arr
Filt = DI.Filt
using .Munge.timeshift: getshift
using .DIutils.statistic: pfunc
filt_desc = Filt.get_filters_desc()

using DataStructures: OrderedDict
import DimensionalData: Between
import DI, DIutils.Table
using ProgressMeter, DimensionalData, Infiltrator, JLD2, DataFrames,
      DataFramesMeta, StatsBase, HypothesisTests, Plots, StatsPlots,
      Statistics, NaNStatistics

import DataStructures: OrderedDict
import DimensionalData: Between
