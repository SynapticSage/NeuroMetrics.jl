using GoalFetchAnalysis
GFA =GoalFetchAnalysis
using GFA.Timeshift, GFA.Plot, GFA.Timeshift.types, GFA.Timeshift.shiftmetrics,
    GFA.Field.metrics, GFA.Plot.receptivefield, GFA.DIutils.namedtup,
    GFA.Munge.isolated, GFA.Munge.nonlocal, GFA.Munge.spiking, GFA.Plot.lfplot,
    GFA.DIutils.arr
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

using GLM, Lasso, Distributions, ThreadSafeDicts, DrWatson, SoftGlobalScope

import DataStructures: OrderedDict
import DimensionalData: Between

using GLMNet, MultivariateStats, MLJ, ScikitLearn
using MATLAB
using Metrics, .Table, GLMNet
import Dates
mat"dbclear all"
using MLJ
using MLJScikitLearnInterface: ElasticNetCVRegressor

function commit_cycwise_vars()::Bool
    has_df = false
    jldopen(path_iso(opt; append="_cyclewise"), "a"; compress=true) do storage
        for name in ("grd","occ","df", "model_isocount","model_hasiso",
                     "model_spikecount", "ca1cycstat", "pfccycstat", 
                     "cyccellstat", "model_cellhasiso", "model_cellhasiso_matlab")
            if isdefined(Main, Symbol(name))
                @info "saving $name"
                obj = @eval Main eval(Symbol($name))
                if name in keys(storage)
                    delete!(storage, name)
                end
                storage[name] = obj
            end
        end
        has_df = "df" in keys(storage)
    end
    has_df
end

# Create a little save function I can run at any time
# (these vars are refs to the data, and if those change,
#  will commit them)
function commit_vars()
    varnames = (
        ("lfp", "spikes", "tsk", "cells", "beh", "cycles", "Rdf", "opt")
    )
    jldopen(path_iso(opt), "a") do storage
        for n in varnames
            v = @eval Main eval(Symbol($n))
            if n in keys(storage)
                delete!(storage, n)
            end
            storage[n] = v
        end
    end
end
