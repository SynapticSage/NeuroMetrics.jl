module fixed

    using DataStructures
    using DataFrames
    import Base
    using LoopVectorization
    using Infiltrator
    using ProgressMeter
    using Entropies: Probabilities
    using ProgressLogging
    using RecipesBase
    using Statistics

    using ..Field
    import ..Field: Grid, ReceptiveField, Occupancy
    import ..Field: get_boundary, resolution_to_width, return_vals
    import ..Field.metrics: MetricSet
    import Utils
    import Table
    import Table: CItype, CItype_plusNull
    import Utils.filtreg: filterAndRegister

    
                                                   

    struct GridFixed <: Field.Grid
        props::Array{String}
        centers::Tuple
        edges::Tuple
        grid::Array{Array{Float32}}
        function GridFixed(props::Vector, centers::Union{Array,Tuple}) 
            if eltype(props) == Symbol
                props = String.(props)
            end
            centers = centers isa Array ? Tuple(centers) : centers
            GridFixed(props,centers)
        end
        function GridFixed(props::Vector, centers::Union{<:AbstractArray,Tuple}, radii::Real)
            if eltype(props) == Symbol
                props = String.(props)
            end
            centers = centers isa Array ? Tuple(centers) : centers
            grid = Field.sparse_to_full(centers)
            radii = ones(size(grid))*radii
            edges = Field.center_to_edge.([[c...] for c in centers])
            edges = Tuple((e...,) for e in edges)
            new(props,centers, edges, grid)
        end
        function GridFixed(props::Vector, centers::Tuple, grid::Array, radii::Array)
            if eltype(props) == Symbol
                props = String.(props)
            end
            @assert(size(grid) == size(radii))
            edges = Field.center_to_edge.([[c...] for c in centers])
            edges = Tuple((e...,) for e in edges)
            new(props,centers,edges,grid)
        end
        function GridFixed(props::Vector;width::Vector, boundary::Vector)
            centers = Tuple(Tuple(collect(s:w:e))
                            for (w, (s, e)) in zip(width, boundary))
            GridFixed(props,centers)
        end
    end

    struct FixedOcc <: Field.Occupancy
        grid::GridFixed
        count::Array{Int32}
        prob::Probabilities
        camerarate::Float32
    end

    struct FixedRF <: ReceptiveField
        grid::GridFixed
        occ::FixedOcc
        count::Array{Int32}
        rate::Array{Float32}
        metrics::MetricSet
    end

    # Setup iteration
    Base.length(g::GridFixed)  = length(g.grid)
    Base.size(g::GridFixed)    = size(g.grid)
    Base.iterate(g::GridFixed) = Base.iterate(g.grid)
    #Base.done(g::GridFixed, state::Int) = length(g.centers) == state
    function Base.iterate(g::GridFixed, state::Tuple{Int,Int})
        iterate(g.grid, state)
    end
    cenumerate(g::GridFixed) = zip(CartesianIndices(g.grid), g)

    FixedFieldDict    = OrderedDict{<:NamedTuple, FixedRF}
    #FixedFieldDict(x) = OrderedDict{NamedTuple,   FixedRF}(x)
    #FixedFieldDict()  = OrderedDict{NamedTuple,   FixedRF}()
    

    """
        find_grid

    obtains the dynamic sampling grid from only the animals behavior
    """
    function get_grid_bounded(behavior, props; 
            #eliminateunusedrowscols::Bool=true,
            widths::OrderedDict, 
            boundary::OrderedDict)::GridFixed
        cv(x) = collect(values(x))
        G = GridFixed(props;width=cv(widths), boundary=cv(boundary))
        G
    end

    function get_grid(behavior::DataFrame, props::Vector;
            widths::Union{<:Int, Vector{<:Int}, OrderedDict},
            other_kws=(;))::GridFixed
        if typeof(widths) <: Int
            widths = OrderedDict{}(prop=>widths for prop in props)
        elseif typeof(widths) <: AbstractVector
            widths = OrderedDict{}(prop=>width for (width,prop) in zip(widths,props))
        else
            @assert(widths isa OrderedDict)
        end
        boundary = get_boundary(behavior, props)
        find_grid_bounded(behavior, props; widths=widths, boundary, other_kws...)
    end
    

    #################################################
    ##### classic fixed grid based  #################
    #################################################

    """
        RFs(spikes, behavior, props; kws...)

    computes an fixed grid and ratemap based on methods in ulanovsky papers
    """
    function RFs(spikes::DataFrame, behavior::DataFrame, props::Vector;
            splitby::Vector=nothing, grid_kws...)::Union{FixedFieldDict, 
                                                 FixedRF}
        grid = find_grid(behavior, props; grid_kws...)
        occ  = get_occupancy(behavior, props, grid)
        RFs(spikes, grid, occ; splitby, grid_kws...)
    end
    function RFs(spikes::DataFrame, grid::GridFixed, occ::FixedOcc;
            splitby::CItype_plusNull=[:unit],
            grid_kws...)::Union{FixedFieldDict, FixedRF}
        if splitby !== nothing
            spikes = groupby(spikes, splitby)
            fields = get_adaptivefields(spikes, grid, occ)
        else
            fields = get_adaptivefield(spikes, grid, occ)
        end
        return fields 
    end

    """
        get_adaptivefields(spikeGroups::GroupedDataFrame, props::Vector,
        grid::GridFixed; kws...)::FixedFieldDict

    computes fixed ratema based on a fixed grid derived from behavior
    """
    function get_fixedfields(spikeGroups::GroupedDataFrame, 
            grid::GridFixed, occ::FixedOcc)::FixedFieldDict
        D = OrderedDict{NamedTuple, FixedRF}()
        #V = []
        keys_and_groups = collect(zip(Table.group.nt_keys(spikeGroups),
                                      spikeGroups))
        #Prog = Progress(length(keys_and_groups); desc="units")
        for (nt, group) in keys_and_groups
            D[nt] = get_adaptivefield(DataFrame(group), grid, occ)
        end
        #for (i, (nt, group)) in enumerate(keys_and_groups)
        #    @infiltrate
        #    D[nt] = fetch(V[i])
        #end
        return D
    end

    """
        get_adaptivefield(X::DataFrame, props::Vector,
                          grid::GridFixed, occ::FixedOcc)::FixedRF

    computes fixed ratemap based on a fixed grid derived from behavior
    """
    function get_fixedfield(spikes::DataFrame, 
            grid::GridFixed, occ::FixedOcc)::FixedRF
        vals = return_vals(spikes, props)
        count = zeros(Int32, size(grid))
        #prog = Progress(length(grid))
        Threads.@threads for (index, (center, radius)) in collect(enumerate(grid))
            count[index] = sum(inside(vals, center, radius))
            #next!(prog)
        end
        count = reshape(count, size(grid))
        FixedRF(grid, occ, count, occ.camerarate*Float32.(count./occ.count))
    end

    function get_occupancy(behavior::DataFrame, 
            grid::GridFixed)::FixedOcc
        vals = return_vals(behavior, props)
        count = zeros(Int32, size(grid))
        prog = Progress(length(grid);desc="occupancy")
        Threads.@threads for (index, (center, radius)) in collect(enumerate(grid))
            count[index] = sum(inside(vals, center, radius))
            next!(prog)
        end
        count = reshape(count, size(grid))
        prob = Probabilities(Float32.(vec(count)))
        camerarate = Float32(1/median(diff(behavior.time)))
        if camerarate > 300
            camerarate /= 60
        end
        FixedOcc(grid, count, prob, camerarate)
    end

    # Helper fucntions
    function vector_dist(vals::Array, center::Array)::Vector{Float32}
        sqrt.(sum((vals .- center[Utils.na, :]).^2,dims=2)[:,1])
    end
    function inside(vals::Array, center::Array, radius::Float32)::BitVector
        vector_dist(vals, center) .< radius
    end
    function get_samptime(vals::Array, center::Array, radius::Float32;
            sampletime::Float32=1/30)::Float32
        sum(inside(vals, center, radius)) * sampletime
    end

    ## --------
    ## UTILITIES
    ## --------
    function to_dict(F::FixedRF)
        FF = Dict{Symbol, Any}
        FF[:count]        = F.count
        FF[:rate]         = F.rate
        FF[:occ_prob]     = reshape(F.occ.prob, size(F.grid.grid))
        FF[:occ_count]    = F.occ.count
        FF[:grid_centers] = F.grid.centers
        FF[:props]        = F.grid.props
        FF
    end

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
    #    G = GridFixed(width, boundary, width)
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

