module adaptive

    using  ..Field
    import ..Field: Grid, ReceptiveField, Occupancy
    import ..Field: get_boundary, resolution_to_width, return_vals
    import ..Field.metrics: Metrics, push_metric!, pop_metric!
    import ..Field: metrics
    import Utils
    import Table
    import Filt

    using DataStructures
    using DataFrames
    import Load.utils: filterAndRegister, register
    import Base
    using LoopVectorization
    using ProgressLogging, ProgressMeter
    using Entropies: Probabilities
    using ProgressLogging
    using RecipesBase
    using Statistics
    using Polyester
    using ThreadSafeDicts
    using Infiltrator

    metric_def = [metrics.bitsperspike, metrics.totalcount, metrics.maxrate,
                  metrics.maxcount, metrics.meanrate, metrics.coherence,
                  metrics.centroid, metrics.argmax, metrics.convexhull]
    
    using Plots
    using DataFrames: ColumnIndex
    CItype = Union{ColumnIndex, Vector{<:ColumnIndex}}
    CItype_plusNull = Union{ColumnIndex, Vector{<:ColumnIndex}, Nothing}

    export yartsev

    """
        default_radii

    Determines how much to expand a hyperphere to fit the corners of the grid
    tagent to the sphere
    """
    function default_radii(centers::Tuple)::Union{Float32, Vector{Float32}}
        C = Vector{Float32}()
        # Push the max diff of each dimension
        for center in centers
            c = maximum(diff([center...]))
            push!(C, c)
        end
        # Turn diameter into radius
        C ./= 2
        # Return the euclidean distance to a corner of the hypercube from the
        # center
        if length(unique(C)) == 1
            sqrt(sum(C.^2))
        else
            C
        end
    end

    struct GridAdaptive <: Grid

        props::Array{String}
        centers::Tuple
        edges::Tuple
        grid::Array{Array{Float32}}
        radii::Union{Array{Float32},
                     Array{Vector{Float32}}}

        function GridAdaptive(props::Vector, centers::Union{Array,Tuple}) 
            if eltype(props) == Symbol
                props = String.(props)
            end
            centers = centers isa Array ? Tuple(centers) : centers
            radii = default_radii(centers)
            GridAdaptive(props,centers, radii)
        end
        function GridAdaptive(props::Vector, centers::Union{<:AbstractArray,Tuple}, 
                radii::Union{Float32, Vector{Float32}})
            if eltype(props) == Symbol
                props = String.(props)
            end
            centers = centers isa Array ? Tuple(centers) : centers
            centers = Tuple(Float32.(c) for c in centers)
            grid = Field.sparse_to_full(centers)
            radii = fill(radii, size(grid))
            edges = Field.center_to_edge.([[c...] for c in centers])
            edges = Tuple((e...,) for e in edges)
            new(props, centers, edges, grid, radii)
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
            GridAdaptive(props, centers)
        end
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

    
    """
        converge_to_radius

    simple converge to radius. it expands a circle/sphere until a certain
    amount of sample time is acheived. it stops when threshold is reached.
    if it reaches maxrad radius before that threshold, this radius becomes
    NaN, nullifying this sample.
    """
    function converge_to_radius(vals::Matrix{Float32}, center::Vector{Float32},
            radius::Float32; dt::Float32, thresh::Float32,
            maxrad::Float32, radiusinc::Float32, kws...)::Float32
        #list = []
        #push!(list, radius)
        @inbounds while get_samptime(vals, center, radius; dt) < thresh
            @fastmath radius = radius * radiusinc
            #push!(list, radius)
            if radius > maxrad
                radius = NaN32
                break
            end
        end
        #@infiltrate
        radius
    end
    function converge_to_radius(vals::Matrix{Float32}, center::Vector{Float32},
            radius::Vector{Float32}; dt::Float32, thresh::Float32,
            maxrad::Union{Float32,Vector{Float32}},
            radiusinc::Union{Float32,Vector{Float32}}, kws...)::Vector{Float32}
        #list = []
        #push!(list, radius)
        #radius = copy(radius)
        @inbounds while get_samptime(vals, center, radius; dt) < thresh
            #@info radius
            @fastmath radius = radius .* radiusinc
            #push!(list, radius)
            if any(radius .> maxrad) || any(radius .=== NaN32)
                radius .= NaN32
                break
            end
        end
        #@infiltrate
        radius
    end

    function converge_to_radius_w_inertia(vals::Matrix{Float32},
            center::Vector{Float32}, radius::Float32; dt::Float32,
            thresh::Float32, maxrad::Float32, radiusinc::Float32,
            ϵ::Float32, kws...)::Float32
        Δ = []
        tₛ = get_samptime(vals, center, radius; dt)
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
                tₛ = get_samptime(vals, center, radius; dt)
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
        get_grid_bounded

    obtains the dynamic sampling grid from only the animals behavior

    ### Notes
    if radiusinc is vector, then it will instead expand a hyperbox instead 
    of a hypersphere
    """
    function get_grid_bounded(behavior::DataFrame, props::Vector; 
            thresh::Float32=1.25f0, # Threshold in seconds
            dt::Union{Nothing,Float32}=nothing, # Total time of sample
            radiusinc::Union{Float32,Vector{Float32}}=0.2f0, # Spatial unit of RF
            ϵ::Float32=0.1f0,
            maxrad::Union{Float32,Nothing}=nothing,
            method::Symbol=:converge_to_radius,
            info::Bool=false,
            widths::OrderedDict, boundary::OrderedDict, kws...)::GridAdaptive
        behavior = Field.ensureTimescale(behavior)
        dt = dt === nothing ? 
                     median(diff(behavior.time)) : dt
        widths = OrderedDict(prop=>widths[prop] for prop in props)
        maxrad = maxrad === nothing ? 3*maximum(values(widths)) : maxrad
        vals = Field.return_vals(behavior, props)
        cv(x) = collect(values(x))
        G = GridAdaptive(props; width=cv(widths), boundary=cv(boundary))
        radiusinc = ((valtype(G.radii) <: Vector) && !(typeof(radiusinc) <: Vector)) ?
                fill(radiusinc, length(G.centers)) : radiusinc
        R = Vector{valtype(G.radii)}(undef, length(G))
        if method == :converge_to_radius_w_inertia
            thresh += ϵ
        end
        method = eval(method)
        radiusinc = radiusinc .+ 1
        if info || !(isdefined(Main, :PlutoRunner))
            @info "grid" thresh dt maxrad radiusinc 
            @info widths
            prog = P = Progress(length(G); desc="Grid")
        end
        #i = 0
        Threads.@threads for (index, (center, radius)) in collect(enumerate(G))
            #i+=1
            #@info "pre $(i)" radius G.radii
            newradius = method(vals, center, radius; dt, thresh, maxrad,
                            radiusinc, ϵ)
            #@info "post $(i)" newradius G.radii
            #@infiltrate
            R[index] = newradius
            if !(isdefined(Main, :PlutoRunner))
                next!(P)
            end
        end
        G.radii .= reshape(R, size(G))
        return G
    end

    function get_grid(behavior::DataFrame, props::Vector;
            widths::Union{<:Float32, Vector{<:Float32}, OrderedDict},
            other_kws...)::GridAdaptive
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
    thread_field_default  = true
    thread_fields_default = false

    """
        yartsev(spikes, behavior, props; kws...)

    computes an adaptive grid and ratemap based on methods in yartsev papers
    """
    function yartsev(behavior::DataFrame, spikes::DataFrame, props::Vector;
            splitby::CItype_plusNull=[:unit], 
            filters::Union{<:AbstractDict, Nothing}=nothing, 
            metrics::Union{Function, Vector{Function}, Nothing}=metric_def, 
            thread_field::Bool=thread_field_default,
            thread_fields::Bool=thread_fields_default,
            grid_kws...)::Union{AdaptiveFieldDict, AdaptiveRF}
        if filters !== nothing
            if Filt.filters_use_precache(filters) &&
                Filt.missing_precache_output_cols(spikes, filters)
                @info "yartsev preacaching"
                spikes = Filt.precache(spikes, behavior, filters)
            end
            @info "filter and reg"
            @time behavior, spikes = filterAndRegister(behavior, spikes; filters,
                                                 on="time", transfer=props,
                                                 filter_skipmissingcols=true)
        else
            exists = intersect(Symbol.(props), 
                                 propertynames(spikes))
            if length(exists) != length(props)
                behavior, spikes = register(behavior, spikes; 
                                            on="time",transfer=String.(props))
                                                     
            end
        end
        @info "grid"
        @time grid = get_grid(behavior, props; grid_kws...)
        @info "occupancy"
        @time occ  = get_occupancy(behavior, grid)
        @info "dropmissing"
        @time spikes = dropmissing(spikes, props)
        @info "fields"
        @time yartsev(spikes, grid, occ; splitby, metrics, grid_kws...)
    end
    function yartsev(spikes::DataFrame, grid::GridAdaptive, occ::AdaptiveOcc;
            splitby::CItype=[:unit],
            metrics::Union{Function, Vector{Function}, Nothing}=metric_def,
            thread_field::Bool=thread_field_default,
            thread_fields::Bool=thread_fields_default,
            grid_kws...)::Union{AdaptiveFieldDict, AdaptiveRF}

        get_adaptivefields(groupby(spikes, splitby), grid, occ;
                                    metrics, thread_field, thread_fields)
    end

    """
        get_adaptivefields(spikeGroups::GroupedDataFrame, props::Vector,
        grid::GridAdaptive; kws...)::AdaptiveFieldDict

    computes adaptive ratema based on a fixed grid derived from behavior
    """
    function get_adaptivefields(spikeGroups::GroupedDataFrame, 
            grid::GridAdaptive, occ::AdaptiveOcc; 
            thread_fields::Bool=thread_fields_default,
            metrics::Union{Function, Vector{Function}, Nothing}=metric_def,
        kws...)::AdaptiveFieldDict
        keys_and_groups = collect(zip(Table.group.nt_keys(spikeGroups),
                                      spikeGroups))
        if thread_fields
            D = ThreadSafeDict{NamedTuple, AdaptiveRF}()
            Threads.@threads for (nt, group) in keys_and_groups
                D[nt] = get_adaptivefield(DataFrame(group), grid, occ; metrics, kws...)
            end
            D = OrderedDict(D)
        else
            D = OrderedDict{NamedTuple, AdaptiveRF}()
            for (nt, group) in keys_and_groups
                D[nt] = get_adaptivefield(DataFrame(group), grid, occ; metrics, kws...)
            end
        end
        return D
    end

    """
        get_adaptivefield(X::DataFrame, props::Vector,
                          grid::GridAdaptive, occ::AdaptiveOc)::AdaptiveRF

    computes adaptive ratemap based on a fixed grid derived from behavior
    """
    function get_adaptivefield(spikes::DataFrame, 
            grid::GridAdaptive, occ::AdaptiveOcc;
            metrics::Union{Function, Vector{Function}, Nothing}=metric_def,
            thread_field::Bool=thread_field_default,
        )::AdaptiveRF
        vals = Field.return_vals(spikes, grid.props)
        count = zeros(Int32, size(grid))
        if thread_field
            #@info "thread single field"
            Threads.@threads for (index, (center, radius)) in
                collect(enumerate(grid))
                @inbounds count[index] = sum(inside(vals, center, radius))
            end
        else
            @inbounds for (index, (center, radius)) in
                collect(enumerate(grid))
                count[index] = sum(inside(vals, center, radius))
            end
        end
        count = reshape(count, size(grid))
        rate  = @fastmath occ.camerarate*Float32.(count./occ.count)
        field = AdaptiveRF(grid, occ, count, rate, Metrics())
        if metrics !== nothing
            for metric in metrics
                push_metric!(field, metric)
            end
        end
        field
    end

    function get_occupancy(behavior::DataFrame, grid::GridAdaptive)::AdaptiveOcc
        vals = Field.return_vals(behavior, grid.props)
        count = zeros(Int32, size(grid))
        Threads.@threads for (index, (center, radius)) in collect(enumerate(grid))
            @inbounds count[index] = sum(inside(vals, center, radius))
        end
        count = reshape(count, size(grid))
        prob  = Probabilities(Float32.(vec(count)))
        camerarate = Float32(1/median(diff(behavior.time)))
        if camerarate > 300
            camerarate /= 60
        end
        AdaptiveOcc(grid, count, prob, camerarate)
    end

    # ----------------
    # Helper functions
    # ----------------
    # Radius measurements
    function vector_dist(vals::Array, center::Array)::Vector{Float32}
        @inbounds @fastmath sqrt.(sum((vals .- center[Utils.na, :]).^2,dims=2)[:,1]) # HUGE SAVINGS FAST MATH HERE
    end
    function inside(vals::Array, center::Array, radius::Float32)::BitVector
        vector_dist(vals, center) .< radius
    end
    function get_samptime(vals::Array, center::Array, radius::Float32;
            dt::Float32=1/30)::Float32
        @fastmath sum(inside(vals, center, radius)) * dt
    end
    # Vector radius measurments
    function indiv_dist(vals::Array, center::Array)::Matrix{Float32}
        abs.(vals .- center[Utils.na, :])
    end
    function inside(vals::Array, center::Array, radius::Vector{Float32})::BitVector
        Utils.squeeze(all(indiv_dist(vals, center) .< radius[Utils.na, :], dims=2))
    end
    function get_samptime(vals::Array, center::Array, radius::Vector{Float32};
            dt::Float32=1/30)::Float32
        @fastmath sum(inside(vals, center, radius)) * dt
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
    #        dt=1/30, # Total time of sample
    #        radiusinc=0.1, # Spatial unit of RF
    #        width::OrderedDict, boundary::OrderedDict)
    #    vals_behavior = return_vals(behavior, props)
    #    vals_spikes   = return_vals(spikes, props)
    #    G = GridAdaptive(width, boundary, width)
    #    @avx for (index, center, radius) in cenumerate(G)
    #        while (sum(vals .< (center .+ radius)) * dt) < thresh
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
