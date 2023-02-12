module Munge

    using Revise
    using Reexport
    using DrWatson
    import ..GoalFetchAnalysis: Timeshift

    include("Munge/chrono.jl")
    include("Munge/behavior.jl")
    include("Munge/lfp.jl")
    include("Munge/well.jl")
    include("Munge/task.jl")
    include("Munge/tensor.jl")
    include("Munge/spiking.jl")
    include("Munge/timeshift.jl")
    include("Munge/lfp_decode.jl")
    #includer"Munge/SpikeTrains.jl")
    #includer"Munge/fieldgrad.jl")
    include("Munge/manifold.jl")
    include("Munge/causal.jl")
    include("Munge/triggering.jl")
    include("Munge/isolated.jl")
    #include(srcdir("Munge","dynamic.jl"))
    include("Munge/nonlocal.jl")

end

