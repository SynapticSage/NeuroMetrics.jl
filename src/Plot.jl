module Plot

    using DrWatson
    using Reexport
    using Plots
    using Infiltrator

    folder_args = []
    exts = ["png", "pdf"]

    function setfolder(args::String...)
        [pop!(folder_args) for arg in 1:length(folder_args)]
        [push!(folder_args, item) for item in args]
        folder = plotsdir(folder_args...)
        !(isdir(folder)) ? mkpath(folder) : nothing
        folder
    end
    function path_to_folderargs(arg::String)
        standard = plotsdir()
        @assert occursin(standard, arg)
        args = replace(arg, standard=>"")
        if startswith(args, "/")
            args[2:end]
        else
            args
        end
    end
    setfolder_from_path(arg::String) = setfolder(path_to_folderargs(arg))

    function save(desc::String; rmexist=nothing)
        folder = plotsdir(folder_args...)
        core = joinpath(folder, desc)
        names = [join([core, ext], ".") for ext in exts]
        for name in names
            @info "saving" name
            savefig(name)
        end
        Plots.CURRENT_PLOT
    end
    function save(desc::NamedTuple; linker="=>", delim=",", rmexist=nothing)
        desc = ["$(k)$(linker)$(v)" for (k,v) in zip(keys(desc),values(desc))]
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

    include(srcdir("Plot","raster.jl"))
    @reexport using .raster

    include(srcdir("Plot","timeshift.jl"))
    @reexport using .timeshift

    include(srcdir("Plot","receptivefield.jl"))
    @reexport using .receptivefield

    include(srcdir("Plot","table.jl"))
    @reexport using .table

    include(srcdir("Plot","notebook_compareTS.jl"))
    @reexport using .notebook_compareTS

    export raster, table

end
