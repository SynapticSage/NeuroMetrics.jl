# General Utilities
using DrWatson
using LinearAlgebra
using ProgressMeter
using Dates
using Logging
using MATLAB

# Data Manipulation and Statistics
using DataFrames, DataFramesMeta, StatsBase, HypothesisTests, Statistics, NaNStatistics
using Distributions, GLM, Lasso, GLMNet, MultivariateStats
using DimensionalData: Between
using DataStructures: OrderedDict
using ThreadSafeDicts
using DIutils.Table

# Plotting and Visualization
using Plots, StatsPlots

# Distributed Computing
using Distributed

# Project-Specific Modules
using GoalFetchAnalysis
import GoalFetchAnalysis.Field.metrics: skipnan
using DI, DIutils.namedtup, DIutils.arr
using DIutils.statistic: pfunc
using .Timeshift, .Timeshift.types, .Timeshift.shiftmetrics
using .Field.metrics
using .Munge.isolated, .Munge.nonlocal, .Munge.spiking
import .Munge.timeshift: getshift
using .Plot.lfplot
using .Plot, .Plot.receptivefield

# Miscellaneous
using SoftGlobalScope

mat"dbclear all"
# using MLJ
# using MLJScikitLearnInterface: ElasticNetCVRegressor # BUG: temporarily causes segfault

export commit_vars
"""
    commit_vars()
Save the variables in the global scope to a jld2 file. The variables must be in
the global scope, and must be one of the following:
lfp, spikes, tsk, cells, beh, cycles, Rdf, opt
"""
function commit_vars(vars::Union{Nothing,Vector,Tuple,String}=nothing)
    if "commits" ∉ keys(opt); opt["commits"] = false; end
    if !opt["commits"]; 
        println("Not committing variables to file.")
        return Nothing; end
    varnames = (x for x in
        ("lfp", "spikes", "tsk", "cells", "beh", "cycles", "Rdf", "opt","R")
        if isdefined(Main, Symbol(x))
    )
    varnames = vars === nothing ? varnames : intersect(varnames, vars)
    jldopen(GoalFetchAnalysis.Munge.isolated.path_iso(opt), "a") do storage
        for n in varnames
            v = @eval Main eval(Symbol($n))
            if n in keys(storage)
                delete!(storage, n)
            end
            storage[n] = v
        end
    end
end

"""
    print_vars()

Print the variables in path_iso(opt; append="")
"""
function print_vars()
    jldopen(path_iso(opt; append=""), "r") do storage
        println(keys(storage))
    end
end

export commit_cycwise_vars
"""
    commit_cyclewise_vars(vars::Union{Nothing,Vector,Tuple,String}=nothing)

Save the variables in the global scope to a jld2 file.

# Arguments
- vars: a vector of strings, or a string, or nothing. If nothing, save all
variables in the global scope. If a string, save that variable. If a vector
of strings, save those variables.
The variables must be in the global scope, and must be one of the following:
grd, occ, df, ca1cycstat, pfccycstat, cyccellstat, model_*, shuffle_*,
glm_list*, glm_dict*
"""
function commit_cycwise_vars(vars::Union{Nothing,Vector,Tuple,String}=nothing)
    has_df = false
    if !in("commits", keys(opt)); opt["commits"] = true; end
    if opt["commits"] == false; 
        println("Not committing variables to file.")
        return Nothing; end
    jldopen(path_iso(opt; append="_cyclewise"), "a"; compress=true) do storage
        models = [string(x) for x in propertynames(Main) 
                  if occursin("model_", string(x))]
        shuffles = [string(x) for x in propertynames(Main) 
                  if occursin("shuffle_", string(x))]
        result_dumps = [string(x) for x in propertynames(Main) 
                  if occursin("glm_list", string(x)) || 
                     occursin("glm_dict", string(x))]
        total_vars = ("grd","occ","df", "ca1cycstat", "pfccycstat",
            "cyccellstat", models..., shuffles..., result_dumps...)
        vars = vars isa String ? [vars] : vars
        vars = vars===nothing ? total_vars : 
                collect(intersect(total_vars, vars))
        @info "Var list to save $vars"
        @showprogress "saving" for name in vars
            if isdefined(Main, Symbol(name)) && 
                !isa(getproperty(Main, Symbol(name)), Function) &&
                !isa(getproperty(Main, Symbol(name)), Module)
                !isa(getproperty(Main, Symbol(name)), Type)

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

"""
    print_cyclewise_vars()
Print the variables in path_iso(opt; append="_cyclewise")
# Arguments
- li: a string, if provided, only print variables that contain `li`
- showsize: a boolean, if true, print the size of each variable
""" 
function print_cyclewise_vars(li=nothing;showsize::Bool=false)
    li = li isa String ? [li] : li
    jldopen(path_iso(opt; append="_cyclewise"), "r") do storage
        if li === nothing
            K = keys(storage)
        else
            K = filter(x -> string(x) ∈ li, keys(storage))
        end
        if showsize
            # compute size of each variable
            for k in K
                # get cumulative size of deep object
                println(k, " : ", Base.summarysize(storage[k]))
            end
        else
            println(K)
        end
    end
end

export initorget
"""
    initorget(key; append="_cyclewise", obj=ThreadSafeDict())
Get the object `key` from the jld2 file at path_iso(opt; append="_cyclewise").
If the key is not in the file, return `obj`.
Will only load if opt["overwrite"] is false.
# Arguments
 - key: a string, the key to get
 - append: a string, the string to append to the path_iso(opt) to get the path
   to the jld2 file.
 - obj: an object, the object to return if the key is not in the file.
# Returns
 - the object stored at key in the jld2 file, or `obj` if the key is not in the
   file.
"""
function initorget(key;append="_cyclewise", obj=ThreadSafeDict())
    println("Mode : overwrite=", opt["overwrite"])
    out = jldopen(path_iso(opt;append), "r") do storage
        print("Possible keys", keys(storage))
        out = if !opt["overwrite"] && key in keys(storage)
            @info "Loading $key from file"
            storage[key]
        else
            obj
        end
    end
    if out isa AbstractDict
        summary(out)
    end
    out
end

"""
    run_shuffle!(df, pos...; shuffle_models=ThreadSafeDict(), shufcount=10_000,
        scrub=["ypred", "m"], kws...)

Shuffle the dataframe labels and run the model on the shuffled data.
"""
function run_shuffle!(df, pos...; 
    shuffle_models=ThreadSafeDict(), shufcount=10_000, scrub=["ypred", "m"],
    kws...)
    for i in 1:shufcount
        key = hash(i + rand())
        shuffle_models[key] = Dict()
        run_glm!(shuffle_cyclelabels(df),
            pos...; kws..., modelz=shuffle_models[key])
        if !isempty(scrub)
            shuffle_models[key] = 
                isolated.clean_keys(shuffle_models[key], scrub)
            GC.gc(); 
        end
    end
end

"""
    checkbins(var)

Print the number of unique bins in the variable `var`.
"""
function checkbins(var::Symbol)
    println("Unique $var bins: ", length(unique(getproperty(Main,var).bins)))
end

predkey(k::NamedTuple) = (;k..., bin=0)

@info "finsished imports_isolated"



"""
    dataframe(G::Vector)
Convert a vector of `DataFrame`s, `Vector`s, and `Dict`s to a single
`DataFrame`.
"""
function todataframe(G::Vector)
    DF = DataFrame()
    @showprogress for g in G
        df = if g isa DataFrame
            g
        elseif g isa Vector
            dataframe(g)
        elseif g isa Dict
            dictdf=DataFrame()
            # Convert to a dataframe row
            for k in filter(k->g[k] !== nothing, keys(g))
                try
                if g[k] isa Vector
                    dictdf[!,k] = [g[k]]
                elseif g[k] isa AbstractDict
                    for kk in keys(g[k])
                        dictdf[!,string(k)*"_"*string(kk)] = [g[k][kk]]
                    end
                else
                    dictdf[!,k] = [g[k]]
                end
                catch
                    @info "Error in $k $(g[k])"
                end
            end
            dictdf
        elseif g isa Nothing
            DataFrame()
        end
        append!(DF, df, cols=:union)
    end
    DF
end

# Write a function that takes a dict, and converts each value to
# a Vector{Vector{Value}} if it's not a vector. If it is a vector,
# then it should be a vector of vectors.

"""
    convert_to_vv(d::Dict)
Convert a `Dict` to a `Dict` where each value is a `Vector` of `Vector`s.
# Arguments
 - d: a `Dict`
# Returns
- d: a `Dict` where each value is a `Vector` of `Vector`s.
"""
function convert_to_vv(d::Dict)
    # if any values of the dict are a dict, then merge them with d
    for (k,v) in d
        if isa(v, AbstractDict)
            d = copy(d)
            d = merge(d, pop!(d, k))
        end
    end
    # if any values of the dict are not a vector, then convert them to a vector
    for (k,v) in d
        if !isa(v, Vector)
            d[k] = [v]
        elseif !isa(v[1], Vector)
            d[k] = [v]
        end
    end
    d
end

