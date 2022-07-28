module Munge

    using Revise
    using Reexport
    using DrWatson

    include(srcdir("Munge","chrono.jl"))
    @reexport using .chrono
    include(srcdir("Munge","behavior.jl"))
    @reexport using .behavior
    include(srcdir("Munge","lfp.jl"))
    @reexport using .lfp
    include(srcdir("Munge","well.jl"))
    @reexport using .well
    include(srcdir("Munge","task.jl"))
    @reexport using .task
    include(srcdir("Munge","tensor.jl"))
    @reexport using .tensor
    include(srcdir("Munge","spiking.jl"))
    @reexport using .spiking

end
