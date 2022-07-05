module fixed

    using ..Field
    using ..Field.RF
    using ..Field.adaptive: AdaptiveOcc, AdaptiveRF, GridAdaptive
    using ..Field: RF, Grid, Occupancy
    import Utils
    import Table
    using DataStructures
    using DataFrames
    import Base
    using LoopVectorization
    using Infiltrator
    using ProgressMeter
    using Entropies: Probabilities
    using ProtoStructs
    using ProgressLogging
    using RecipesBase
    using Statistics

    
                                                   

    struct GridFixed <: Field.Grid
        props::Array{String}
        centers::Tuple
        edges::Tuple
        grid::Array{Array{Float32}}
        function GridAdaptive(props::Vector, centers::Union{Array,Tuple}) 
            if eltype(props) == Symbol
                props = String.(props)
            end
            centers = centers isa Array ? Tuple(centers) : centers
            GridAdaptive(props,centers)
        end
        function GridAdaptive(props::Vector, centers::Union{<:AbstractArray,Tuple}, radii::Real)
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
        function GridAdaptive(props::Vector, centers::Tuple, grid::Array, radii::Array)
            if eltype(props) == Symbol
                props = String.(props)
            end
            @assert(size(grid) == size(radii))
            edges = Field.center_to_edge.([[c...] for c in centers])
            edges = Tuple((e...,) for e in edges)
            new(props,centers,edges,grid)
        end
        function GridAdaptive(props::Vector;width::Vector, boundary::Vector)
            centers = Tuple(Tuple(collect(s:w:e))
                            for (w, (s, e)) in zip(width, boundary))
            GridAdaptive(props,centers)
        end
    end

    struct FixedOcc <: Field.Occupancy
        grid::GridFixed
        count::Array{Int32}
        prob::Probabilities
        camerarate::Float32
    end

    struct FixedRF <: Field.ReceptiveField
        grid::GridAdaptive
        occ::FixedOcc
        count::Array{Int32}
        rate::Array{Float32}
    end

    # Setup iteration
    Base.length(g::GridAdaptive)  = length(g.grid)
    Base.size(g::GridAdaptive)    = size(g.grid)
    Base.iterate(g::GridAdaptive) = Base.iterate(zip(g.grid, g.radii))
    #Base.done(g::GridAdaptive, state::Int) = length(g.centers) == state
    function Base.iterate(g::GridAdaptive, state::Tuple{Int,Int})
        iterate(zip(g.grid,g.radii), state)
    end
    cenumerate(g::GridAdaptive) = zip(CartesianIndices(g.grid), g)

    AdapativFieldDict    = OrderedDict{<:NamedTuple, AdaptiveRF}
    #AdaptiveFieldDict(x) = OrderedDict{NamedTuple,   AdaptiveRF}(x)
    #AdaptiveFieldDict()  = OrderedDict{NamedTuple,   AdaptiveRF}()
    

    """
        find_grid

    obtains the dynamic sampling grid from only the animals behavior
    """
    function find_grid_bounded(behavior, props; 
            thresh::Float32=1f0, # Threshold in seconds
            sampletime::Float32=1/30f0, # Total time of sample
            radiusinc::Float32=0.5f0, # Spatial unit of RF
            ϵ::Float32=0.1f0,
            maxrad::Float32=5f0,
            method::Symbol=:converge_to_radius,
            #eliminateunusedrowscols::Bool=true,
            widths::OrderedDict, boundary::OrderedDict)::GridAdaptive
        vals = return_vals(behavior, props)
        cv(x) = collect(values(x))
        G = GridAdaptive(props;width=cv(widths), boundary=cv(boundary))
        P = Progress(length(G), desc="grid")
        R = Vector{Float32}(undef, length(G))
        if method == :converge_to_radius_w_inertia
            sampletime += ϵ
        end
        method = eval(method)
        Threads.@threads for (index, (center, radius)) in collect(enumerate(G))
            #@infiltrate
            radius = method(vals, center, radius; sampletime, thresh, maxrad,
                            radiusinc, ϵ)
            next!(P)
            R[index] = radius
        end
        G.radii .= reshape(R, size(G))
        #if eliminateunusedrowscols
        #    rows = findall(all(isnan.(G.radii), dims=2))
        #    cols = findall(all(isnan.(G.radii), dims=1))
        #    G.grid = G.grid[Not(rows), Not(cols)]
        #end
        G
    end

    function find_grid(behavior::DataFrame, props::Vector;
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
        find_grid_bounded(behavior, props; widths=widths, boundary, other_kws...)
    end
    

    #################################################
    ##### ulanovsky paper based  ####################
    #################################################

    """
        ulanovsky(spikes, behavior, props; kws...)

    computes an adaptive grid and ratemap based on methods in ulanovsky papers
    """
    function RFs(spikes::DataFrame, behavior::DataFrame, props::Vector;
            splitby::Vector=nothing, grid_kws...)::Union{AdapativFieldDict, 
                                                 AdaptiveRF}
        grid = find_grid(behavior, props; grid_kws...)
        occ  = get_occupancy(behavior, props, grid)
        ulanovsky(spikes, props, grid, occ; splitby, grid_kws...)
    end
    function RFs(spikes::DataFrame, props::Vector, grid::GridAdaptive,
            occ::AdaptiveOcc;
            splitby::Vector=nothing, grid_kws...)::Union{AdapativFieldDict, 
                                                 AdaptiveRF}
        if splitby !== nothing
            spikes = groupby(spikes, splitby)
            fields = get_adaptivefields(spikes, props, grid, occ)
        else
            fields = get_adaptivefield(spikes, props, grid, occ)
        end
        return fields 
    end

    """
        get_adaptivefields(spikeGroups::GroupedDataFrame, props::Vector,
        grid::GridAdaptive; kws...)::AdapativFieldDict

    computes adaptive ratema based on a fixed grid derived from behavior
    """
    function get_fixedfields(spikeGroups::GroupedDataFrame, props::Vector,
            grid::GridAdaptive, occ::AdaptiveOcc)::AdapativFieldDict
        D = OrderedDict{NamedTuple, AdaptiveRF}()
        #V = []
        keys_and_groups = collect(zip(Table.group.nt_keys(spikeGroups),
                                      spikeGroups))
        #Prog = Progress(length(keys_and_groups); desc="units")
        for (nt, group) in keys_and_groups
            D[nt] = get_adaptivefield(DataFrame(group), props, grid, occ)
        end
        #for (i, (nt, group)) in enumerate(keys_and_groups)
        #    @infiltrate
        #    D[nt] = fetch(V[i])
        #end
        return D
    end

    """
        get_adaptivefield(X::DataFrame, props::Vector,
                          grid::GridAdaptive, occ::AdaptiveOcc)::AdaptiveRF

    computes adaptive ratemap based on a fixed grid derived from behavior
    """
    function get_fixedfield(spikes::DataFrame, props::Vector,
            grid::GridAdaptive, occ::AdaptiveOcc)::AdaptiveRF
        vals = return_vals(spikes, props)
        count = zeros(Int32, size(grid))
        #prog = Progress(length(grid))
        Threads.@threads for (index, (center, radius)) in collect(enumerate(grid))
            count[index] = sum(inside(vals, center, radius))
            #next!(prog)
        end
        count = reshape(count, size(grid))
        AdaptiveRF(grid, occ, count, occ.camerarate*Float32.(count./occ.count))
    end

    function get_occupancy(behavior::DataFrame, props::Vector,
            grid::GridAdaptive)::AdaptiveOcc
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
        AdaptiveOcc(grid, count, prob, camerarate)
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

    ## --------
    ## UTILITIES
    ## --------
    function to_dict(F::FixedRF)
        FF = Dict{Symbol, Any}
        FF[:count] = F.count
        FF[:rate] = F.rate
        FF[:occ_prob] = F.occ.prob
        FF[:occ_count] = F.occ.count
        FF[:grid_centers] = F.grid.centers
        FF[:grid] = F.grid.grid
        FF
    end

end

