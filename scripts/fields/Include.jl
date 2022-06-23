using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))

# Grab our raw data
using DataFrames
using KernelDensity, Distributions
using Plots, Measures
using ProgressMeter
using StatsPlots
using DataFramesMeta
using DataStructures: OrderedDict
using Revise
using Distributions

using GoalFetchAnalysis
import Field
import Load
#utils = GoalFetchAnalysis.utils ALREADY EXPORTED IN USING ABOVE

operation     = Field.operation
model         = Field.model
recon_process = Field.recon_process
