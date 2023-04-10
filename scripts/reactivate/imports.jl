
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
using Arrow
using HypothesisTests
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
    vars = vars isa String ? [vars] : vars
    if !in("commits", keys(opt)); opt["commits"] = true; end
    if opt["commits"] == false; 
        println("Not committing variables to file.")
    return Nothing; end
    if vars === nothing || "DF" in vars && isdefined(Main, :DF)
        @info "Committing $(path_react(opt))"
        arrow_file = replace(path_react(opt),"jld2"=>"arrow")
        Arrow.write(arrow_file, DF, compress=:lz4)
    end
    if vars === nothing || "DFS" in vars && isdefined(Main, :DFS)
        @info "Committing $(path_react(opt;append="_summary"))"
        arrow_file = replace(path_react(opt;append="_summary"),"jld2"=>"arrow")
        Arrow.write(arrow_file, DFS, compress=:lz4)
    end
end

def = isdefined(Main, :opt)
opt = def ? opt : Dict()
begin
    if !def
        opt["animal"] = "RY16"
        opt["day"] = 36
    end
end

