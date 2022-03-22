module GV

    using DrWatson
    using Plots

    include(srcdir("raw.jl"))
    include(srcdir("table.jl"))
    include(srcdir("raster.jl"))
    include(srcdir("goalvector.jl"))
    export raw 
    export table 
    export raster 
    export goalvector
    export Plots

end
