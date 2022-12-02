module task

    using RecipesBase
    using DataFrames
    using Infiltrator
    using LazyGrids
    using LazySets
    using Colors

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

    @userplot Heatmap_NoBoundary
    @recipe function heatmap_noboundary(plt::Heatmap_NoBoundary; 
                                        resolution=80,
                                        epoch=nothing, transpose=false)
        task     = plt.args[1]
        epoch    = (epoch === nothing) ? 
                    first(task[task.type .== "run",:].epoch) : epoch
        boundary = task[(task.name.=="boundary") .& 
                        (task.epoch .== epoch), :]
        append!(boundary, DataFrame(boundary[1,:]))

        label     --> "boundary"
        c         --> :orange
        linestyle --> :dash
        #seriestype := :heatmap
        
        # Create a grid
        x = LinRange(minimum(task.x), maximum(task.x), resolution)
        y = LinRange(minimum(task.y), maximum(task.y), resolution)
        X,Y = ndgrid(x,y);
        #@infiltrate

        # Create a hull with the boundary
        H = VPolygon(hcat(boundary.x,boundary.y)';apply_convex_hull=false)
        #plot(eachrow(hcat(H.vertices...))...)

        # Z = (X,Y) ∈ Hull
        Z = [element(s) ∈ H for s in Singleton.(collect.(zip(X,Y)))]
        @assert any(Z)
        Zim = [notin==1 ? RGBA(1,1,1,1) : RGBA(NaN,NaN,NaN,0) for notin in Z]
        seriesalpha --> (!).(Z)
        #seriestype := :image
        
        !transpose ? (x,y,Zim') : (y,x,Zim)
    end

end

