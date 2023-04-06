
# This script grabs reactivation scores with my Munge.reactivation module
using GoalFetchAnalysis
using GoalFetchAnalysis.Munge.reactivation, GoalFetchAnalysis.Munge.spiking,
      GoalFetchAnalysis.Munge.groupop
using DI
using DimensionalData
using DIutils.filtreg
import DIutils
using DataFrames, ProgressMeter, DataFramesMeta, StatsPlots
using GoalFetchAnalysis.Munge.spiking, DIutils.arr, Infiltrator, Statistics
using Plots, DataVoyager, Logging, JLD2
import DataStructures: OrderedDict
LinearAlgebra.BLAS.set_num_threads(16)
using MKL
ENV["MKL_NUM_THREADS"] = 16
export commit_react_vars
"""
    commit_react_vars(vars::Union{Nothing,Vector,Tuple,String}=nothing)
Save the variables in the global scope to a jld2 file.
# Arguments
- vars: a vector of strings, or a string, or nothing. If nothing, save all
variables in the global scope. If a string, save that variable. If a vector
of strings, save those variables.
The variables must be in the global scope, and must be one of the following:
grd, occ, df, ca1cycstat, pfccycstat, cyccellstat, model_*, shuffle_*,
glm_list*, glm_dict*
"""
function commit_react_vars(vars::Union{Nothing,Vector,Tuple,String}=nothing)
    has_df = false
    if !in("commits", keys(opt)); opt["commits"] = true; end
    if opt["commits"] == false; 
        println("Not committing variables to file.")
        return Nothing; end
    jldopen(path_react(opt), "a"; compress=true) do storage
        total_vars = ("Scores", "DF", "props_ca1", "props_pfc")
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
    end
end


def = isdefined(Main, :opt)
begin
    opt = def ? opt : Dict()
    if !def
        opt["animal"] = "RY16"
        opt["day"] = 36
    end
end

