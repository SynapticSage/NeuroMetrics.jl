module Plot

    using DrWatson
    using Reexport
    using Plots
    using Infiltrator

    folder_args = []
    exts = ["png", "pdf"]

    function save(desc::String)
        core = plotsdir(folder_args..., desc)
        names = [join([core, ext], ".") for ext in exts]
        for name in names
            @info "saving" name
            savefig(name)
        end
        Plots.CURRENT_PLOT
    end
    function save(desc::NamedTuple; linker="=>", delim=",")
        desc = ["$(k)$(linker)$(v)" for (k,v) in zip(keys(desc),values(desc))]
        desc = join(desc, delim)
        save(desc)
    end
    function save(desc::Dict; kws...)
        save(NamedTuple(desc); kws...)
    end

    include(srcdir("Plot","raster.jl"))
    @reexport using .raster

    include(srcdir("Plot","timeshift.jl"))
    @reexport using .timeshift

    include(srcdir("Plot","receptivefield.jl"))
    @reexport using .receptivefield

    include(srcdir("Plot","table.jl"))
    @reexport using .table

    export raster, table

end
