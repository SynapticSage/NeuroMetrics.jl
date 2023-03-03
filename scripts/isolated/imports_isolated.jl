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

using GLM, Lasso, Distributions, ThreadSafeDicts, DrWatson, SoftGlobalScope

import DataStructures: OrderedDict
import DimensionalData: Between

using GLMNet, MultivariateStats, MLJ, ScikitLearn
using MATLAB
using Metrics, .Table, GLMNet
import Dates, Logging
mat"dbclear all"
using MLJ
using MLJScikitLearnInterface: ElasticNetCVRegressor

export commit_vars
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

export commit_cycwise_vars
"""
    commit_cyclewise_vars

"""
function commit_cycwise_vars()::Bool
    has_df = false
    jldopen(path_iso(opt; append="_cyclewise"), "a"; compress=true) do storage
        models = [string(x) for x in propertynames(Main) 
                  if occursin("model_", string(x))]
        shuffles = [string(x) for x in propertynames(Main) 
                  if occursin("shuffle_", string(x))]
        for name in ("grd","occ","df", "ca1cycstat", "pfccycstat", 
                     "cyccellstat", models..., shuffles...)
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

export initorget
function initorget(key;append="_cyclewise", obj=ThreadSafeDict())
    jldopen(path_iso(opt;append), "r") do storage
        print("Possible keys", keys(storage))
        if key in keys(storage)
            storage[key]
        else
            obj
        end
    end
end


using Random
export shuffle_cyclelabels
function shuffle_cyclelabels(df, perm=[:cycs,:relcycs])
    # if "orig"*string(perm[1]) ∉ names(df)
    #         df[!,"orig".*string(perm)] .= df[!, perm]
    # end
    inds = collect(1:size(df,1))
    randperm!(inds)
    sort(
    DataFrames.transform(df, perm .=> x->x[inds], Not(perm), renamecols=false),
        perm)
end

function run_shuffle!(df, pos...; 
    shuffle_models=ThreadSafeDict(), shufcount=10_000, scrub_ypred=true, 
    scrub_model=true, kws...)
    @showprogress "shuffle" for i in 1:shufcount
        key = hash(i+rand())
        shuffle_models[key] = Dict()
        run_glm!(shuffle_cyclelabels(df),
            pos...; kws..., modelz=shuffle_models[key])
        if scrub_ypred
            [pop!(v, "ypred") for v in values(shuffle_models[key])
            if !(v isa Exception)]
        end
        if scrub_model
            [pop!(v, "m") for v in values(shuffle_models[key])
            if !(v isa Exception)]
        end
        if scrub_ypred || scrub_model; GC.gc(); end
    end
end

function get_dx_dy(df, relcyc)
    dx = @subset(df, :relcycs .== relcyc)
    dy = @subset(df, :relcycs .== 0)
    register = [:cyc_batch, :cyc_match]
    # dx = groupby(dx, [:cyc_batch, :cyc_match])
    # dy = groupby(dy, [:cyc_batch, :cyc_match])
    # kx, ky = keys(dx.keymap), keys(dy.keymap)
    # k = intersect(kx,ky)
    # dx = sort(vcat([dx[kk] for kk in k]...), [:cyc_batch, :cyc_match])
    # dy = sort(vcat([dy[kk] for kk in k]...), [:cyc_batch, :cyc_match])
    dx, dy = sort(dx, register), sort(dy, register)
    dxr, dyr = Matrix(dx[!,register]), Matrix(dy[!,register])
    if dxr != dyr
        dxr, dyr = eachrow.((dxr, dyr))
        D = intersect(dxr, dyr)
        idx = [findfirst((x == d for x in dxr)) for d in D]
        idy = [findfirst((y == d for y in dyr)) for d in D]
        dx, dy = dx[idx,:], dy[idy, :]
    end
    @assert Matrix(dx[!,register]) == Matrix(dy[!,register])
    dx, dy
end

function get_futurepast_blocks(df)
    dxf = unstack(@subset(df, :relcycs .> 0), )
    dxp = @subset(df, :relcycs .<= 0)
    dy = @subset(df, :relcycs .== 0)
    dxp = groupby(dxp, [:cyc_batch, :cyc_match])
    dxf = groupby(dxf, [:cyc_batch, :cyc_match])
    dy = groupby(dy, [:cyc_batch, :cyc_match])
    kxp, kxf, ky = keys(dxp.keymap), keys(dxf.keymap), keys(dy.keymap)
    k = intersect(kxf,kxp,ky)
    dxf = sort(vcat([dxf[kk] for kk in k]...), [:cyc_batch, :cyc_match], )
    dxp = sort(vcat([dxp[kk] for kk in k]...), [:cyc_batch, :cyc_match], )
    dy  = sort(vcat([dy[kk] for kk in k]...), [:cyc_batch, :cyc_match], )
    Dict("pastcurr"=>(dxp, dy), "future"=>(dxf, dy))
end
