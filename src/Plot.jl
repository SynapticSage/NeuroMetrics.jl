module Plot

using DrWatson
using Reexport
using Plots
using Infiltrator
using RecipesBase
using Measures
using Utils
using Infiltrator

parent_folder = []
folder_args = []
exts = ["png", "pdf"]
append, prepend = "", ""

setappend(val)  = @eval Plot append = $val
setprepend(val) = @eval Plot prepend = $val

function setparentfolder(args::String...)
    [pop!(parent_folder) for arg in 1:length(parent_folder)]
    [push!(parent_folder, item) for item in args]
    folder = plotsdir(parent_folder..., folder_args...)
    !(isdir(folder)) ? mkpath(folder) : nothing
    folder
end
function setfolder(args::String...)
    [pop!(folder_args) for arg in 1:length(folder_args)]
    [push!(folder_args, item) for item in args]
    folder = plotsdir(parent_folder..., folder_args...)
    !(isdir(folder)) ? mkpath(folder) : nothing
    folder
end
function path_to_folderargs(arg::String)
    standard = plotsdir()
    @assert occursin(standard, arg)
    args = replace(arg, standard => "")
    if startswith(args, "/")
        args[2:end]
    else
        args
    end
end
setparentfolder_from_path(arg::String) = setparentfolder(path_to_folderargs(arg))
setfolder_from_path(arg::String) = setfolder(path_to_folderargs(arg))

function save(desc::String; rmexist=nothing)
    desc = replace(desc,
                   "OrderedDict" => "",
                   "Dict" => "",
                   "String" => "",
                   "Float32" => "", "Float64" => "", "{," => "", 
                   "{"=>"", "}" => "")
    folder = plotsdir(parent_folder..., folder_args...)
    append_string = append isa NamedTuple ? 
             Utils.namedtup.ntopt_string(append) : append
    prepend_string = prepend isa NamedTuple ? 
             Utils.namedtup.ntopt_string(prepend) : prepend
    names = [join([prepend_string, desc, append_string, ext], ".") for ext in exts]
    names = [startswith(name,".") ? name[2:end] : name  for name in names]
    names = replace.(names, ".."=>".")
    names = joinpath.([folder], names)
    for name in names
        @info "saving" name
        savefig(name)
    end
    Plots.CURRENT_PLOT
end
function save(desc::NamedTuple; linker="=>", delim=",", rmexist=nothing)
    desc = ["$(k)$(linker)$(v)" for (k, v) in zip(keys(desc), values(desc))]
    desc = join(desc, delim)
    save(desc; rmexist)
end
function save(desc::Dict; kws...)
    save(NamedTuple(desc); kws...)
end
function save(plot::Plots.Plot, pos...; kws...)
    save(pos...; kws...)
    plot
end
function save(folder::String, plot::Plots.Plot, pos...; kws...)
    setfolder_from_path(folder)
    save(pos...; kws...)
    plot
end

function create_blank_plot()
    plot(legend=false, grid=false, framestyle=:none, background_color_inside=:match)
end


# -----------------
# Steroscopic plots
# -----------------
export stereoplot
@userplot StereoPlot
@recipe function stereoplot(plt::StereoPlot; theta=0, phi=35, offset=6)
    x,y,z = plt.args 
    seriestype --> :scatter
    margin --> -6mm
    layout := (1,2)
    @series begin
        projection_type := :perspective
        subplot := 1
        camera := (mod(theta,360+offset),phi)
        (x,y,z)
    end
    projection_type := :perspective
    subplot := 2
    camera := (mod(theta,360),phi)
    (x,y,z)
end

function stereoscopicgif(pos...;delta_angle=1,kws...)
    @gif for theta in 1:delta_angle:360
        stereoplot(pos...; theta, kws...)
    end
end

#include(srcdir("Plot", "raster.jl"))
#@reexport using .raster

include(srcdir("Plot", "timeshift.jl"))
@reexport using .timeshift

include(srcdir("Plot", "receptivefield.jl"))
@reexport using .receptivefield

include(srcdir("Plot", "table.jl"))
@reexport using .table

#include(srcdir("Plot", "notebook_compareTS.jl"))
#@reexport using .notebook_compareTS

include(srcdir("Plot", "task.jl"))
@reexport using .task

include(srcdir("Plot", "cause.jl"))
@reexport using .cause

include(srcdir("Plot", "nonlocal.jl"))
@reexport using .nonlocal

end
