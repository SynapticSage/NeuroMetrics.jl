"""
Module for triggered measurements

This module allows for taking samples when data crosses into a radii of measurement

"""
module triggering

    using DataFrames
    using Utils.binning
    using Table: CItype
    using Infiltrator
    import Utils

    WindowType = Union{Real, Tuple{<:Real,<:Real}}
    """
        TriggerGenerator

    Object that instantiates grabbing a window around a trigger
    """
    mutable struct TriggerGenerator
        X::DataFrame # Data triggering on
        grid::GridAdaptive
        occ::IndexedAdaptiveOcc
        win::WindowType
    end
    Base.length(g::TriggerGenerator) = length(g.grid)
    Base.size(g::TriggerGenerator)   = size(g.grid)
    """
    iterators that loops through the trigger generator and outputs
    (grid, data) points
    """
    function Base.getindex(g::TriggerGenerator, state)
        window_centroids = g.X[g.occ.inds[state],"time"]
        windows_of_data = Vector{DataFrame}(undef, size(window_centroids,1))
        for trig in eachrow(window_centroids)
            times = Utils.in_range(g.X.time, g.win .* (-1,1) .+ trig)
            push!(windows_of_data, g.X[times, :])
        end
        (g.grid.grid[state], windows_of_data)
    end
    Base.iterate(g::TriggerGenerator) = g[1], 1
    Base.iterate(g::TriggerGenerator, state) = state >= length(g) ? 
                                        nothing : (g[state+1], state+1)
    Base.peek(g::TriggerGenerator)  = g[1]
    Base.peek(g::TriggerGenerator, state) = g[state]
    #Base.popfirst!(g::TriggerGenerator)   = nothing

    export get_triggergen
    """
        get_triggergen(X::DataFrame, props::CItype, window::WindowType; grid_kws...)

    Gets a trigger generator object
    """
    function get_triggergen(X::DataFrame, props::CItype, window::WindowType; grid_kws...)
        @debug "get_triggergen :: getting grid"
        grid = get_grid(X, props ;grid_kws...)
        @debug "get_triggergen :: getting occ"
        occ  = get_occupancy_indexed(X, grid)
        #fields = union(string.(props), ["time"])
        TriggerGenerator(X, grid, occ, window) 
    end
    function get_triggergen(X::DataFrame, props::CItype, window::WindowType,
            grid::Grid)
        occ = get_occupancy_indexed(X, grid)
        TriggerGenerator(X, grid, occ, win) 
    end
    function get_triggergen(X::DataFrame, props::CItype, window::WindowType,
            grid::Grid, occ::Occupancy)
        TriggerGenerator(X, grid, occ, win) 
    end
end
