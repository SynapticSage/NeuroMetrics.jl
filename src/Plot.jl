module Plot

using DrWatson, Reexport, Plots, Infiltrator, RecipesBase, Measures,
      Infiltrator
using DIutils

parent_folder = []
folder_args   = []
complete_folder_args = []
exts = ["png", "pdf"]
append, prepend = "", ""
active = true

setappend(val)  = @eval Plot append = $val
appendtoappend(val)  = @eval Plot append = append * $val
setprepend(val) = @eval Plot prepend = $val
prependtoprepend(val)  = @eval Plot prepend = $val * prepend
function off()
    @eval Plot active = false
end
function on()
    @eval Plot active = true
end
toggle() = @eval Plot !active

function setparentfolder(args::String...)
    @eval Plot parent_folder = $args
    @eval Plot complete_folder_args = [parent_folder...,folder_args...]

    folder = plotsdir(parent_folder..., folder_args...)
    !(isdir(folder)) ? mkpath(folder) : nothing
    folder
end
function setfolder(args::String...)
    @eval Plot folder_args = $args
    @eval Plot complete_folder_args = [parent_folder...,folder_args...]

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
             DIutils.namedtup.ntopt_string(append) : append
    prepend_string = prepend isa NamedTuple ? 
             DIutils.namedtup.ntopt_string(prepend) : prepend
    names = [join([prepend_string, desc, append_string, ext], ".") for ext in exts]
    names = [startswith(name,".") ? name[2:end] : name  for name in names]
    names = replace.(names, ".."=>".")
    names = joinpath.([folder], names)
    for name in names
        active ? begin 
            @info "saving" name
            savefig(name)
            end : nothing
    end
    Plots.CURRENT_PLOT
end
function save(desc::NamedTuple; linker="=>", delim=",", rmexist=nothing,append=nothing,prepend=nothing)
    desc = ["$(k)$(linker)$(v)" for (k, v) in zip(keys(desc), values(desc))]
    desc = join(desc, delim)
    desc = append === nothing ? desc  : desc * append
    desc = prepend === nothing ? desc : prepend * desc 
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

 include("Plot/timeshift.jl")
 include("Plot/receptivefield.jl")
 include("Plot/table.jl")
 include("Plot/task.jl")
 include("Plot/cause.jl")
 include("Plot/nonlocal.jl")
 include("Plot/lfplot.jl")
 include("Plot/manifold.jl")
 export timeshift, receptivefield, table, task, 
        cause, nonlocal, lfplot, manifold

# include(srcdir("Plot", "timeshift.jl"))
# @reexport using .timeshift
# include(srcdir("Plot", "receptivefield.jl"))
# @reexport using .receptivefield
# include(srcdir("Plot", "table.jl"))
# @reexport using .table
# #include(srcdir("Plot", "notebook_compareTS.jl"))
# #@reexport using .notebook_compareTS
# include(srcdir("Plot", "task.jl"))
# @reexport using .task
# include(srcdir("Plot", "cause.jl"))
# @reexport using .cause
# include(srcdir("Plot", "nonlocal.jl"))
# @reexport using .nonlocal
# include(srcdir("Plot", "lfplot.jl"))
# @reexport using .lfplot
# include(srcdir("Plot", "manifold.jl"))
# @reexport using .manifold

end
