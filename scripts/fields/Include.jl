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
includet(srcdir("raw.jl"))
includet(srcdir("field.jl"))
includet(srcdir("field/operation.jl"))
includet(srcdir("field/model.jl"))
includet(srcdir("field/recon_process.jl"))
includet(srcdir("filt.jl"))
includet(srcdir("table.jl"))
includet(srcdir("utils.jl"))
includet(srcdir("shuffle.jl"))
