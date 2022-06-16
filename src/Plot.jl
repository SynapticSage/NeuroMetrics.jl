module Plot

    include(srcdir("Plot","raster.jl"))
    @reexport using .raster
    include(srcdir("Plot","timeshift.jl"))
    @reexport using .timeshift
    include(srcdir("Plot","table.jl"))
    @reexport using .table

    export raster, timeshift, table

end
