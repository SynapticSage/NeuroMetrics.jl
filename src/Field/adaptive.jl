module adaptive

    using DataStructures
    using DataFrames
    using LazyGrids: ndgrid
    import Base
    using ..Field.RF
    using LoopVectorization

    abstract type Grid end
    abstract type RF end
    abstract type FieldDict end

    grid_to_full(G::Vector...) = [g for g in zip(ndgrid(G...))]
    grid_to_full(G::Tuple)     = [g for g in zip(ndgrid(G...))]
                                                   

    struct GridAdaptive <: Grid
        centers::Array
        radii::Array
        GridAdaptive(c::Array, r::Array) = new(c,r)
        GridAdaptive(c::Tuple, r::Tuple) = new(grid_to_full(c),grid_to_full(r))
        GridAdaptive(c::Tuple) = new(grid_to_full(c), undef)
    end

    # Setup iteration
    Base.iterate(g::GridAdaptive) = Base.iterate(zip(G.centers,G.radii))
    #Base.done(g::GridAdaptive, state::Int) = length(g.centers) == state
    Base.iterate(g::GridAdaptive, state::Int) = Base.iterate(zip(G.centers, G.radii), 
                                                       state)
    cenumerate(G::GridAdaptive) = zip(CartesianIndices(G.centers), G)

    struct AdaptiveField
        grid::GridAdaptive
        field::Array
    end
    const AdapativFieldDict = OrderedDict{<:NamedTuple, <:Any}

    # Utility functions
    """
        resolution_to_width

    converts resolution to width, ie the number of points spread over a
    range of points into the inter-sample width
    """
    function resolution_to_width(resolution::OrderedDict, boundary::OrderedDict)
        width = OrderedDict{String,Any}()
        for prop in props
            width[prop] = boundary[prop] / resolution[prop]
        end
        width
    end

    """
        get_boundary

    gets the boundary of each property
    """
    function get_boundary(behavior::DataFrame, props::Vector)
        boundary   = OrderedDict{String,Any}()
        for props in props
            boundary[prop]   = extrema(Utils.skipnan(behavior[!,prop]))
        end
    end

    """
        return_vals

    return a value matrix/dataframe for the requested
    properties in the DataFrame X
    """
    function return_vals(X::DataFrame, props::Vector)
        vals = hcat([X[!,prop] for prop in props]...)
    end

    ###########################################
    ##### ULANOVSKY BASED  ####################
    ###########################################

    #function ulanovsky_find_grid(behavior, props; width::Int, kws...)::GridAdaptive
    #    width = OrderedDict(prop=>width for prop in props)
    #    ulanovsky_find_grid(behavior, props; width, boundary, kws...)::GridAdaptive
    #end

    #function ulanovsky_find_grid(behavior, props; widths::Vector{<:Int}, kws...)::GridAdaptive
    #    width = OrderedDict(prop=>width for (prop,width) in zip(props,widths))
    #    ulanovsky_find_grid(behavior, props; width, boundary, kws...)
    #end

    #function ulanovsky_find_grid(behavior, props; width::OrderedDict, kws...)::GridAdaptive
    #    boundary = get_boundary(behavior, props)
    #    ulanovsky_find_grid(behavior, props; width, boundary, kws...)
    #end

    #function ulanovsky_find_grid(behavior, props; 
    #        thresh::Real=1, # Threshold in seconds
    #        sampletime=1/30, # Total time of sample
    #        radiusinc=0.1, # Spatial unit of RF
    #        width::OrderedDict, boundary::OrderedDict, kws...)::GridAdaptive
    #    vals = return_vals(behavior, props)
    #    G = GridAdaptive(width, boundary, width)
    #    for (index, center, radius) in cenumerate(G)
    #        while (sum(vals .< (center .+ radius)) * sampletime) < thresh
    #            radius += radiusinc
    #        end
    #        G.radii[index] = radius
    #    end
    #    G
    #end

    #"""
    #    ulanovsky(spikes, props; grid::GridAdaptive)

    #computes adaptive ratemap based on methods in ulanovsky papers
    #"""
    #function ulanovsky(spikes, props; grid::GridAdaptive)
    #    vals = return_vals(spikes, props)
    #    #vals = replace(vals, NaN=>-Inf)
    #    results = zeros(size(GridAdaptive.centers))
    #    indices = CartesianIndices(results)
    #    @avx for (index, center, radius) in zip(indicies, grid.centers, grid.radii)
    #        results[index] = sum((vals .- center) .< radius)
    #    end
    #end

    #"""
    #    ulanovsky(spikes, behavior, props; kws...)

    #computes an adaptive grid and ratemap based on methods in ulanovsky papers
    #"""
    #function ulanovsky(spikes::DataFrame, behavior::DataFrame, props::Vector; splitby=nothing, kws...)::Union{AdapativFieldDict, AdaptiveField}
    #    G = ulanovsky_find_grid(behavior; kws...)
    #    if splitby != nothing
    #        spikes = groupby(spikes, splitby)
    #    end
    #    ulanovsky(spikes, props; G, splitby)
    #end

    #function ulanovsky(spikeGroups::GroupedDataFrame, props::Vector; G::GridAdaptive, kws...)::AdapativFieldDict
    #    D = AdapativFieldDict()
    #    for (nt, group) in zip(Table.group.nt_keys(spikeGroups), spikeGroups)
    #        D[nt] = ulanovsky(group, props; G, splitby)
    #    end
    #end

    ## ------
    ## Skaggs
    ## ------

    #function skaggs_find_grid(behavior, props; width::Int, kws...)
    #    width = OrderedDict(prop=>width for prop in props)
    #    skaggs_find_grid(behavior, props; width, boundary, kws...)
    #end
    #function skaggs_find_grid(behavior, props; widths::Vector{<:Int}, kws...)
    #    width = OrderedDict(prop=>width for (prop,width) in zip(props,widths))
    #    skaggs_find_grid(behavior, props; width, boundary, kws...)
    #end
    #function skaggs_find_grid(behavior, props; width::OrderedDict, kws...)
    #    boundary = get_boundary(behavior, props)
    #    skaggs_find_grid(behavior, props; width, boundary, kws...)
    #end
    #function skaggs_find_grid(spikes, behavior, props; 
    #        thresh::Real=1, # Threshold in seconds
    #        sampletime=1/30, # Total time of sample
    #        radiusinc=0.1, # Spatial unit of RF
    #        width::OrderedDict, boundary::OrderedDict)
    #    vals_behavior = return_vals(behavior, props)
    #    vals_spikes   = return_vals(spikes, props)
    #    G = GridAdaptive(width, boundary, width)
    #    @avx for (index, center, radius) in cenumerate(G)
    #        while (sum(vals .< (center .+ radius)) * sampletime) < thresh
    #            radius += radiusinc
    #        end
    #        G.radii[index] = radius
    #    end
    #    G
    #end


    #"""
    #"""
    #function skaggs()
    #end

end
