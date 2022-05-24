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
filt = GoalFetchAnalysis.filt
utils = GoalFetchAnalysis.utils

includet(srcdir("raw.jl"))
includet(srcdir("field.jl"))
includet(srcdir("filt.jl"))

operation = field.operation
model = field.model
recon_process = field.recon_process
