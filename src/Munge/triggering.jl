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
        X::Vector{DataFrame} # Data triggering on
        grid::GridAdaptive
        occ::IndexedAdaptiveOcc
        win::T where T<:WindowType
    end
    Base.length(g::TriggerGenerator) = length(g.grid)
    Base.size(g::TriggerGenerator)   = size(g.grid)
    #function TriggerGenerator(X::Vector, grid::GridAdaptive, occ::IndexedAdaptiveOcc, win::T where T<:WindowType)
    #    TriggerGenerator(X,grid,occ,win)
    #end
    """
    iterators that loops through the trigger generator and outputs
    (grid, data) points
    """
    function Base.getindex(g::TriggerGenerator, state)
        window_centroids = g.X[g.occ.inds[state],"time"]
        windows_of_data = Vector{DataFrame}(undef, size(window_centroids,1))
        for trig in eachrow(window_centroids)
            times = Utils.in_range(g.X.time, g.win .* (-1,1) .+ trig)
            push!(windows_of_data, [x[times, :] for x in g.X])
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
    function get_triggergen(props::CItype, window::WindowType, X::DataFrame...; grid_kws...)
        @debug "get_triggergen :: getting grid"
        grid = get_grid(X[1], props ;grid_kws...)
        @debug "get_triggergen :: getting occ"
        occ  = get_occupancy_indexed(X[1], grid)
        #fields = union(string.(props), ["time"])
        TriggerGenerator(collect(X), grid, occ, window) 
    end
    function get_triggergen(props::CItype, window::WindowType, grid::Grid, X::DataFrame...)
        occ = get_occupancy_indexed(X[1], grid)
        TriggerGenerator(collect(X), grid, occ, win) 
    end
    function get_triggergen(props::CItype, window::WindowType, grid::Grid, occ::Occupancy, X::DataFrame...)
        TriggerGenerator(collect(X), grid, occ, win) 
    end

end
