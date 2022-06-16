module Munge
    using Revise
    using Reexport
    using DrWatson

    include(srcdir("Munge","behavior.jl"))
    @reexport using .behavior
    include(srcdir("Munge","lfp.jl"))
    @reexport using .lfp
    include(srcdir("Munge","raster.jl"))
    @reexport using .raster
    include(srcdir("Munge","well.jl"))
    @reexport using .well
end
