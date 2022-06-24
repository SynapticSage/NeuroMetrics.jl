module Plot

    using DrWatson
    using Reexport

    include(srcdir("Plot","raster.jl"))
    @reexport using .raster
    include(srcdir("Plot","timeshift.jl"))
    @reexport using .timeshift
    include(srcdir("Plot","table.jl"))
    @reexport using .table

    function hello3()
        "🍎"
    end

    export raster, table

end
