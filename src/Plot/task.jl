module task

    using RecipesBase
    using DataFrames
    using Infiltrator

    export PlotBoundary
    @userplot PlotBoundary
    @recipe function plot_boundary(plt::PlotBoundary; epoch=nothing, transpose=false)
        task = plt.args[1]
        #@assert (task isa DataFrame) && any("boundary" .== task.name)
        epoch = (epoch === nothing) ? first(task[task.type .== "run",:].epoch) : epoch
        boundary = task[(task.name.=="boundary") .& (task.epoch .== epoch), :]
        append!(boundary, DataFrame(boundary[1,:]))
        label --> "boundary"
        c --> :orange
        linestyle --> :dash
        transpose ? (boundary.x, boundary.y) : (boundary.y, boundary.x)
    end

end

