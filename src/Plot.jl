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
leave_plot_on = true

setappend(val::String)  = (@eval Plot append = $val; ps())
setappend(val::NamedTuple)  = setappend(DIutils.namedtup.ntopt_string(val))
function setappend(val::AbstractDict)
    setappend(DIutils.namedtup.ntopt_string(NamedTuple(Symbol(k)=>v 
    for (k,v) in val)))
end
appendtoappend(val)  = @eval Plot append = append * $val
setprepend(val::String) = (@eval Plot prepend = $val; ps())
setprepend(val::NamedTuple)  = setprepend(DIutils.namedtup.ntopt_string(val))
prependtoprepend(val)  = @eval Plot prepend = $val * prepend
function off()
    @eval Plot active = false
    ps()
end
function on()
    @eval Plot active = true
    ps()
end
toggle() = (@eval Plot !active; ps())

function setparentfolder(args::String...)
    @eval Plot parent_folder = $args
    @eval Plot complete_folder_args = [parent_folder...,folder_args...]

    folder = DrWatson.plotsdir(parent_folder..., folder_args...)
    !(isdir(folder)) ? mkpath(folder) : nothing
    ps()
    folder
end
function setfolder(args::String...)
    @eval Plot folder_args = $args
    @eval Plot complete_folder_args = [parent_folder...,folder_args...]

    folder = DrWatson.plotsdir(parent_folder..., folder_args...)
    !(isdir(folder)) ? mkpath(folder) : nothing
    ps()
    folder
end
function pushfolder(args::String...)
    @eval Plot folder_args = [folder_args..., $args...]
    @eval Plot complete_folder_args = [parent_folder...,folder_args...]

    folder = DrWatson.plotsdir(parent_folder..., folder_args...)
    !(isdir(folder)) ? mkpath(folder) : nothing
    ps()
    folder
end
function popfolder(n::Int=1)
    @eval Plot folder_args = folder_args[1:end-$n]
    @eval Plot complete_folder_args = [parent_folder...,folder_args...]
    folder = DrWatson.plotsdir(parent_folder..., folder_args...)
    !(isdir(folder)) ? mkpath(folder) : nothing
    ps()
    folder
end
function path_to_folderargs(arg::String)
    standard = DrWatson.plotsdir()
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

function inspect(args::String...)
    folder=DrWatson.plotsdir(Plot.parent_folder..., Plot.folder_args..., 
        args...)
    run(`ranger $folder`)
end

function save(desc::String; rmexist=nothing)
    desc = replace(desc,
                   "OrderedDict" => "",
                   "Dict" => "",
                   "String" => "",
                   "Float32" => "", "Float64" => "", "{," => "", 
                   "{"=>"", "}" => "")
    folder = DrWatson.plotsdir(parent_folder..., folder_args...)
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
                if leave_plot_on
                    current()
                end
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

function create_blank_plot(pos...;kws...)
    plot(pos...;legend=false, grid=false, framestyle=:none, background_color_inside=:match, kws...)
end
blank(pos...;kws...) = create_blank_plot(pos...;kws...)


"""
    deleteplotfolder()

Deletes the current plot folder
"""
function deleteplotfolder()
    folder = DrWatson.plotsdir(parent_folder..., folder_args...)
    rm(folder, recursive=true)
    mkpath(folder)
end

"""
    deleteplotfiles()

Deletes all files in the current plot folder
"""
function deleteplotfiles()
    folder = DrWatson.plotsdir(parent_folder..., folder_args...)
    # Prompt user if corret folder
    println("Are you sure you want to delete all files in $folder? (y/n)")
    if readline() != "y"
        return
    else
        @info "Deleting all files in $folder" readdir(folder)
    end
    rm.(joinpath.(folder, readdir(folder)), force=true)
end

function reportplotsinfolder()
    folder = DrWatson.plotsdir(parent_folder..., folder_args...)
    println("Plots in folder: $folder")
    files = readdir(folder)
    # Constrain files to png, pdf, svg
    files = filter(x -> occursin(r"\.(png|pdf|svg)$", x), files)
    println(files)
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

function plotzdir()
    DrWatson.plotsdir(complete_folder_args...)
end

"""
    printstate()

Prints the current state of the plot module, the key variables

- If plots are active
- The plotdir pointed to by the module
- The current append and prepend strings
"""
function printstate()
    println("Plot module state:")
    println("------------------")
    println("Active: ", active)
    println("Plotdir: ", Plot.plotzdir())
    println("Append: ", append)
    println("Prepend: ", prepend)
    println("------------------")
end
ps = printstate

"""
    savetopowerpointslide(file::String, 
                    slide::Int, desc...; rmexist=false)

Saves the current plot to a powerpoint slide.
"""
function savetopowerpointslide(file::String, slide::Int, desc...; 
    rmexist=false)
    folder = DrWatson.plotsdir(parent_folder..., folder_args...)
    append_string = append isa NamedTuple ? 
             DIutils.namedtup.ntopt_string(append) : append
    prepend_string = prepend isa NamedTuple ? 
             DIutils.namedtup.ntopt_string(prepend) : prepend
    name = join([prepend_string, file, append_string, "pptx"], ".")
    name = startswith(name,".") ? name[2:end] : name
    name = replace(name, ".."=>".")
    name = joinpath(folder, name)
    if rmexist === true
        rm(name)
    end
    active ? begin 
        @info "saving" name
        # TODO https://asml-labs.github.io/PPTX.jl/dev/
        end : nothing
    Plots.CURRENT_PLOT
end

function suptitle(suptitle, P=current())
    # Create a subplot with the title but no axes
    suptitle_subplot = Plots.plot(legend=false, grid=false)
    Plots.plot!(title=suptitle, framestyle=:none)

    if P isa Vector
        P = Plots.plot(P...)
    end

    # Combine the suptitle plot with the original plot
    P_new = Plots.plot(suptitle_subplot, P, 
        layout=@layout([a{0.025h}; b{0.975h}]), 
        margin=0Plots.mm)

    return P_new
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
