module adaptive

    using ..Field
    import ..Field: Grid, ReceptiveField, Occupancy, Metrics
    import ..Field: get_boundary, resolution_to_width, return_vals
    import Utils
    import Table
    import Table: CItype, CItype_plusNull

    using DataStructures
    using DataFrames
    import Load.utils: filterAndRegister
    import Base
    using LoopVectorization
    using Infiltrator
    using ProgressMeter
    using Entropies: Probabilities
    using ProgressLogging
    using RecipesBase
    using Statistics
    

    function max_radii(centers::Tuple)
        C = Vector{Float32}()
        for center in centers
            c = maximum(diff([center...]))
            push!(C, c)
        end
        C ./= 2
        sqrt(sum(C.^2))
    end

    struct GridAdaptive <: Grid
        props::Array{String}
        centers::Tuple
        edges::Tuple
        grid::Array{Array{Float32}}
        radii::Array{Float32}
        function GridAdaptive(props::Vector, centers::Union{Array,Tuple}) 
            if eltype(props) == Symbol
                props = String.(props)
            end
            centers = centers isa Array ? Tuple(centers) : centers
            radii = max_radii(centers)
            GridAdaptive(props,centers, radii)
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
            new(props,centers, edges, grid, radii)
        end
        function GridAdaptive(props::Vector, centers::Tuple, grid::Array, radii::Array)
            if eltype(props) == Symbol
                props = String.(props)
            end
            @assert(size(grid) == size(radii))
            edges = Field.center_to_edge.([[c...] for c in centers])
            edges = Tuple((e...,) for e in edges)
            new(props,centers,edges,grid,radii)
        end
        function GridAdaptive(props::Vector; width::Vector, boundary::Vector)
            centers = Tuple(Tuple(collect(s:w:e))
                            for (w, (s, e)) in zip(width, boundary))
            GridAdaptive(props,centers)
        end
    end

    @recipe function plot_adaptiveocc(grid::GridAdaptive, val::Symbol=:radii)
        colorbar_title --> String(val)
        seriestype --> :heatmap
        c --> :thermal
        x --> [grid.centers[1]...]
        if length(grid.centers) > 1
            y --> [grid.centers[2]...]
        end
        getproperty(grid, val)
    end

    struct AdaptiveOcc <: Occupancy
        grid::GridAdaptive
        count::Array{Int32}
        prob::Probabilities
        camerarate::Float32
    end

    struct AdaptiveRF <: ReceptiveField
        grid::GridAdaptive
        occ::AdaptiveOcc
        count::Array{Int32}
        rate::Array{Float32}
        metrics::Metrics
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

    AdaptiveFieldDict    = OrderedDict{<:NamedTuple, AdaptiveRF}
    #AdaptiveFieldDict(x) = OrderedDict{NamedTuple,   AdaptiveRF}(x)
    #AdaptiveFieldDict()  = OrderedDict{NamedTuple,   AdaptiveRF}()
    
    function converge_to_radius(vals::Matrix{Float32}, center::Vector{Float32},
            radius::Float32; sampletime::Float32, thresh::Float32,
            maxrad::Float32, radiusinc::Real, kws...)::Float32
        while get_samptime(vals, center, radius; sampletime) < thresh
            radius += radiusinc
            if radius > maxrad
                radius = NaN
                break
            end
        end
        radius
    end

    function converge_to_radius_w_inertia(vals::Matrix{Float32},
            center::Vector{Float32}, radius::Float32; sampletime::Float32,
            thresh::Float32, maxrad::Float32, radiusinc::Float32,
            ϵ::Float32, kws...)::Float32
        Δ = []
        tₛ = get_samptime(vals, center, radius; sampletime)
        δ = tₛ - thresh
        push!(Δ, δ)
        s = sign(δ)
        tolerance_acheived = abs(δ) < ϵ
        if s == 1
            tolerance_acheived =true
        end
        increase_phase = true
        reversal = 1
        while !(tolerance_acheived)
            while s != reversal && !(tolerance_acheived)
                radius += radiusinc
                if increase_phase 
                    radiusinc *= 2
                else
                    radiusinc /= 2
                end
                tₛ = get_samptime(vals, center, radius; sampletime)
                δ = tₛ - thresh
                push!(Δ, δ)
                s = sign(δ)
                tolerance_acheived = abs(δ) < ϵ
            end
            reversal *= -1
            radiusinc *= -1
            increase_phase = false
            if abs(radiusinc) < 0.01
                break
            end
        end
        #@info "tolernace acheived = $tolerance_acheived"
        if radius > maxrad
            radius = NaN
        end
        radius
    end

    """
        get_grid

    obtains the dynamic sampling grid from only the animals behavior
    """
    function get_grid_bounded(behavior, props; 
            thresh::Float32=1.25f0, # Threshold in seconds
            sampletime::Union{Nothing,Float32}=nothing, # Total time of sample
            radiusinc::Float32=0.5f0, # Spatial unit of RF
            ϵ::Float32=0.1f0,
            maxrad::Float32=5f0,
            method::Symbol=:converge_to_radius,
            widths::OrderedDict, boundary::OrderedDict)::GridAdaptive
        sampletime = sampletime === nothing ? 
                     1/median(diff(behavior.time)) : sampletime
        vals = Field.return_vals(behavior, props)
        cv(x) = collect(values(x))
        G = GridAdaptive(props;width=cv(widths), boundary=cv(boundary))
        P = Progress(length(G), desc="grid")
        R = Vector{Float32}(undef, length(G))
        if method == :converge_to_radius_w_inertia
            sampletime += ϵ
        end
        method = eval(method)
        Threads.@threads for (index, (center, radius)) in collect(enumerate(G))
            radius = method(vals, center, radius; sampletime, thresh, maxrad,
                            radiusinc, ϵ)
            next!(P)
            R[index] = radius
        end
        G.radii .= reshape(R, size(G))
        return G
    end

    function get_grid(behavior::DataFrame, props::Vector;
            widths::Union{<:Float32, Vector{<:Float32}, OrderedDict},
            other_kws=(;))::GridAdaptive
        if typeof(widths) <: Float32
            widths = OrderedDict{}(prop=>widths for prop in props)
        elseif typeof(widths) <: AbstractVector
            widths = OrderedDict{}(prop=>width for (width,prop) in zip(widths,props))
        else
            @assert(widths isa OrderedDict)
        end
        boundary = Field.get_boundary(behavior, props)
        get_grid_bounded(behavior, props; widths=widths, boundary, other_kws...)
    end
    

    #################################################
    ##### yartsev paper based  ####################
    #################################################

    """
        yartsev(spikes, behavior, props; kws...)

    computes an adaptive grid and ratemap based on methods in yartsev papers
    """
    function yartsev(behavior::DataFrame, spikes::DataFrame, props::Vector;
            splitby::CItype_plusNull=[:unit], 
            filters::Union{<:AbstractDict, Nothing}=nothing, 
            grid_kws...)::Union{AdaptiveFieldDict, AdaptiveRF}
        if filters !== nothing
            behavior, spikes = filterAndRegister(behavior, spikes; filters,
                                                 on="time",transfer=props,
                                                 filter_skipmissingcols=true)
        end
        grid = get_grid(behavior, props; grid_kws...)
        occ  = get_occupancy(behavior, grid)
        spikes = dropmissing(spikes, props)
        yartsev(spikes, grid, occ; splitby, grid_kws...)
    end
    function yartsev(spikes::DataFrame, grid::GridAdaptive, occ::AdaptiveOcc;
            splitby::CItype_plusNull=[:unit], 
            grid_kws...)::Union{AdaptiveFieldDict, AdaptiveRF}
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
        grid::GridAdaptive; kws...)::AdaptiveFieldDict

    computes adaptive ratema based on a fixed grid derived from behavior
    """
    function get_adaptivefields(spikeGroups::GroupedDataFrame, 
            grid::GridAdaptive, occ::AdaptiveOcc)::AdaptiveFieldDict
        D = OrderedDict{NamedTuple, AdaptiveRF}()
        #V = []
        keys_and_groups = collect(zip(Table.group.nt_keys(spikeGroups),
                                      spikeGroups))
        #Prog = Progress(length(keys_and_groups); desc="units")
        for (nt, group) in keys_and_groups
            #push!(V, Dagger.spawn(get_adaptivefield, DataFrame(group), props, grid, occ))
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
                          grid::GridAdaptive, occ::AdaptiveOcc)::AdaptiveRF

    computes adaptive ratemap based on a fixed grid derived from behavior
    """
    function get_adaptivefield(spikes::DataFrame, 
            grid::GridAdaptive, occ::AdaptiveOcc)::AdaptiveRF
        vals = Field.return_vals(spikes, grid.props)
        count = zeros(Int32, size(grid))
        #prog = Progress(length(grid))
        Threads.@threads for (index, (center, radius)) in collect(enumerate(grid))
            count[index] = sum(inside(vals, center, radius))
            #next!(prog)
        end
        count = reshape(count, size(grid))
        AdaptiveRF(grid, occ, count, occ.camerarate*Float32.(count./occ.count), 
                   Metrics())
    end

    function get_occupancy(behavior::DataFrame, grid::GridAdaptive)::AdaptiveOcc
        vals = Field.return_vals(behavior, grid.props)
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
    ## --------
    ## UTILITIES
    ## --------
    function Utils.to_dict(F::AdaptiveRF)
        FF = Dict{Symbol, Any}()
        FF[:count]        = F.count
        FF[:rate]         = F.rate
        FF[:occ_prob]     = reshape(F.occ.prob, size(F.grid.grid))
        FF[:occ_count]    = F.occ.count
        FF[:grid_centers] = F.grid.centers
        FF[:props]        = F.grid.props
        FF[:radii]        = F.grid.radii
        FF
    end

    ## ------
    ## Skaggs
    ## ------
    #function skaggs_get_grid(behavior, props; width::Int, kws...)
    #    width = OrderedDict(prop=>width for prop in props)
    #    skaggs_get_grid(behavior, props; width, boundary, kws...)
    #end
    #function skaggs_get_grid(behavior, props; widths::Vector{<:Int}, kws...)
    #    width = OrderedDict(prop=>width for (prop,width) in zip(props,widths))
    #    skaggs_get_grid(behavior, props; width, boundary, kws...)
    #end
    #function skaggs_get_grid(behavior, props; width::OrderedDict, kws...)
    #    boundary = get_boundary(behavior, props)
    #    skaggs_get_grid(behavior, props; width, boundary, kws...)
    #end
    #function skaggs_get_grid(spikes, behavior, props; 
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
