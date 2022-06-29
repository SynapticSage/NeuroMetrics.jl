module Munge

    using Revise
    using Reexport
    using DrWatson

    include(srcdir("Munge","behavior.jl"))
    @reexport using .behavior
    include(srcdir("Munge","lfp.jl"))
    @reexport using .lfp
    include(srcdir("Munge","well.jl"))
    @reexport using .well
    include(srcdir("Munge","task.jl"))
    @reexport using .task

end
