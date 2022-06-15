module Munge
    using Revise
    using Reexport
    using DrWatson

    push!(LOAD_PATH,  srcdir("Munge","src"))
    @reexport using behavior
    @reexport using lfp
    @info LOAD_PATH
    @assert isfile(srcdir("Munge","src","well.jl"))
    @reexport using raster
    @info "got here"
    pop!( LOAD_PATH)
end
