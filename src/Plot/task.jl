module task

    using RecipesBase
    using DataFrames
    using Infiltrator

    export PlotBoundary

    @userplot PlotBoundary
    @recipe function plot_boundary(plt::PlotBoundary; epoch=nothing)
        task = plt.args[1]
        #@assert (task isa DataFrame) && any("boundary" .== task.name)
        boundary = task[(task.name.=="boundary") .& (task.epoch .== epoch), :]
        epoch = (epoch === nothing) ? first(task[task.type .== "run",:].epoch) : epoch
        append!(boundary, DataFrame(boundary[1,:]))
        label --> "boundary"
        c --> :orange
        linestyle --> :dash
        @info boundary
        (boundary.x, boundary.y)
    end

end

