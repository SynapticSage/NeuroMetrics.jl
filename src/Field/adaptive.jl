module adaptive

    using DataStructures
    using DataFrames
    using LazyGrids: ndgrid
    import Base
    using ..Field
    using ..Field.RF
    using LoopVectorization
    using Infiltrator
    import Utils
    using ProgressMeter
    #using ThreadsX

    sparse_to_full(C::Vector...) = [g for g in zip(ndgrid(G...)...)]
    sparse_to_full(C::Tuple{Vector}) = [g for g in zip(ndgrid(G...)...)]
    function sparse_to_full(C::Tuple...)::Array{Tuple}
        C = [[c...] for c in C]
        C = [c for c in zip(ndgrid(C...)...)]
        [collect(c) for c in C]
    end
    function sparse_to_full(C::Tuple)::Array{Array}
        C = [[c...] for c in C]
        C = [c for c in zip(ndgrid(C...)...)]
        return [collect(x) for x in C]
    end
    function full_to_sparse(G::Array)::Array{Array}
        accessor = Vector{Union{Colon,<:Int}}(undef,ndims(G))
        accessor .= 1
        out = []
        for i in 1:ndims(G)
            access = copy(accessor)
            access[i] = Colon()
            g = G[access...]
            g = Tuple(g[i] for g in g)
            push!(out,g)
        end
        out
    end

    function max_radii(centers::Tuple)
        C = Vector{Float64}()
        for center in centers
            c = maximum(diff([center...]))
            push!(C, c)
        end
        C ./= 2
        sqrt(sum(C.^2))
    end
                                                   

    struct GridAdaptive <: Field.Grid
        centers::Tuple
        edges::Tuple
        grid::Array
        radii::Array
        samptime::Array
        function GridAdaptive(pos...; width::Vector, boundary::Vector)
            centers = Tuple(Tuple(collect(s:w:e))
                            for (w, (s, e)) in zip(width, boundary))
            GridAdaptive(centers)
        end
        function GridAdaptive(centers::Union{Array,Tuple}, radii::Real)
            centers = centers isa Array ? Tuple(centers) : centers
            grid = sparse_to_full(centers)
            radii = ones(size(grid))*radii
            edges = Field.center_to_edge.([[c...] for c in centers])
            edges = Tuple((e...,) for e in edges)
            new(centers, edges, grid, radii)
        end
        function GridAdaptive(centers::Union{Array,Tuple}) 
            centers = centers isa Array ? Tuple(centers) : centers
            radii = max_radii(centers)
            GridAdaptive(centers, radii)
        end
        function GridAdaptive(centers::Tuple, grid::Array, radii::Array)
            @assert(size(grid) == size(radii))
            edges = Field.center_to_edge.([[c...] for c in centers])
            edges = Tuple((e...,) for e in edges)
            new(centers,edges,grid,radii)
        end
    end
    # TODO make plot recipe

    # Setup iteration
    Base.length(g::GridAdaptive)  = length(g.grid)
    Base.size(g::GridAdaptive)  = size(g.grid)
    Base.iterate(g::GridAdaptive) = Base.iterate(zip(g.grid,g.radii))
    #Base.done(g::GridAdaptive, state::Int) = length(g.centers) == state
    function Base.iterate(g::GridAdaptive, state::Tuple{Int,Int})
        iterate(zip(g.grid,g.radii), state)
    end
    cenumerate(g::GridAdaptive) = zip(CartesianIndices(g.grid), g)

    struct AdaptiveField
        grid::GridAdaptive
        field::Array
    end
    AdapativFieldDict = OrderedDict{<:NamedTuple, <:Any}

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
        for prop in props
            boundary[prop]   = extrema(Utils.skipnan(behavior[!,prop]))
        end
        boundary
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


    """
        ulanovsky_find_grid

    obtains the dynamic sampling grid from only the animals behavior
    """
    function ulanovsky_find_grid_bounded(behavior, props; 
            thresh::Real=1, # Threshold in seconds
            sampletime=1/30, # Total time of sample
            radiusinc=0.1, # Spatial unit of RF
            maxrad=5,
            widths::OrderedDict, boundary::OrderedDict)::GridAdaptive
        vals = return_vals(behavior, props)
        cv(x) = collect(values(x))
        G = GridAdaptive(;width=cv(widths), boundary=cv(boundary))
        function vector_dist(center,radius)
            sqrt.(sum((vals .- center[Utils.na, :]).^2,dims=2)[:,1])
        end
        inside(center,radius) = vector_dist(center, radius) .< radius
        function get_samptime(center, radius)
            sum(inside(center,radius)) * sampletime
        end
        P = Progress(length(G))
        R = Vector{Float32}(undef, length(G))
        Threads.@threads for (index, (center, radius)) in collect(enumerate(G))
            while get_samptime(center, radius) < thresh
                radius += radiusinc
                #@info (;samptime=get_samptime(center,radius), radiuschange=radius)

                if radius > maxrad
                    radius = NaN
                    break
                end
            end
            next!(P)
            #@info get_samptime(center, radius)
            R[index] = radius
        end
        G.radii .= reshape(R, size(G))
        G
    end
    function ulanovsky_find_grid(behavior::DataFrame, props::Vector;
            widths::Union{<:Int, Vector{<:Int}, OrderedDict},
            other_kws=(;))::GridAdaptive
        if typeof(widths) <: Int
            widths = OrderedDict{}(prop=>widths for prop in props)
        elseif typeof(widths) <: AbstractVector
            widths = OrderedDict{}(prop=>width for (width,prop) in zip(widths,props))
        else
            @assert(widths isa OrderedDict)
        end
        boundary = get_boundary(behavior, props)
        ulanovsky_find_grid_bounded(behavior, props; widths=widths, boundary, other_kws...)
    end
    function test_New()
    end

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
