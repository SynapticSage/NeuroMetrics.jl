"""
Module for triggered measurements
"""
module triggering

    WindowType = Union{Real, Tuple{<:Real,<:Real}}
    struct TrigGrid
        props::Array{String}
        grid::Array{Vector}
        radii::Union{Array, Array{Array}, Real}
    end
    mutable struct TriggerGenerator
        grid::TrigGrid 
        X::DataFrame # Data triggering on
        i::Int # Current state
    end

    function collect_triggers(X::DataFrame, window::WindowType, location::Vector, radius)
    end
    function collect_triggers(X::DataFrame, window::WindowType, grid::TrigGrid)
    end
end
