module Plot

    using DrWatson
    using Reexport

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
