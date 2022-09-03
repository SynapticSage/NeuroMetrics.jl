module Munge

    using Revise
    using Reexport
    using DrWatson

    include(srcdir("Munge","chrono.jl"))
    include(srcdir("Munge","behavior.jl"))
    include(srcdir("Munge","lfp.jl"))
    include(srcdir("Munge","well.jl"))
    include(srcdir("Munge","task.jl"))
    include(srcdir("Munge","tensor.jl"))
    include(srcdir("Munge","spiking.jl"))
    include(srcdir("Munge","timeshift.jl"))
    include(srcdir("Munge","lfp_decode.jl"))
    include(srcdir("Munge","nonlocal.jl"))
    @reexport using .chrono
    @reexport using .behavior
    @reexport using .lfp
    @reexport using .well
    @reexport using .task
    @reexport using .tensor
    @reexport using .spiking
    @reexport using .timeshift
    @reexport using .lfp_decode
    @reexport using .nonlocal

end

