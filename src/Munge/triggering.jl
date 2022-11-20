"""
Module for triggered measurements

This module allows for taking samples when data crosses into a radii of measurement

"""
module triggering

    using Utils.binning
    using DataFrames

    WindowType = Union{Real, Tuple{<:Real,<:Real}}
    mutable struct TriggerGenerator
        X::DataFrame # Data triggering on
        grid::GridAdaptive
        occ::IndexedAdaptiveOcc
        win::WindowType
        i::Int # Current state
    end
    Base.length(g::TriggerGenerator) = length(g.grid)
    Base.size(g::TriggerGenerator)   = size(g.grid)
    """
    iterators that loops through the trigger generator and outputs
    (grid, data) points
    """
    function Base.iterate(g::TriggerGenerator)
        DF_chunks = (g.X[occ.inds[i]] for i in 1:length(occ.inds))
        Base.iterate(zip(g.grid.grid,DF_chunks))
    end
    
    

    export get_triggergen
    function get_triggergen(X::DataFrame, props, window::WindowType, grid_kws...)
        grid = get_grid(X, props ;grid_kws...)
        occ  = get_occupancy_indexed(X, grid)
        TriggerGenerator(X, grid, occ, win, i) 
    end
    function get_triggergen(X::DataFrame, props, window::WindowType, grid::Grid)
        occ = get_occupancy_indexed(X, grid)
        TriggerGenerator(X, grid, occ, win, i) 
    end
    function get_triggergen(X::DataFrame, props, window::WindowType, grid::Grid, occ::Occupancy)
        TriggerGenerator(X, grid, occ, win, i) 
    end




end
