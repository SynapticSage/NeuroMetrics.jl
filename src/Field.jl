module Field

    export get_fields
    export ReceptiveField, Grid, Occupancy


    # Julia Packages
    using DataFrames
    using ImageFiltering
    import DataStructures: OrderedDict
    using DataStructures
    using Statistics
    using NaNStatistics
    using Infiltrator
    using GeometricalPredicates: inpolygon
    using LazyGrids: ndgrid
    using LazySets
    using Statistics

    # Goal Vector Libraries
    using DrWatson
    import Utils
    import Table
    import Load

    rateConversion = 30
    export rateConversion


    # Our basic RECEPTIVE-FIELD type components
    abstract type Grid end
    abstract type ReceptiveField end
    abstract type Occupancy end

    function Base.push!(R::ReceptiveField, measurement_name::Symbol, measurement_value)
        push!(R.metrics, measurement_name => measurement_value)
    end
    function Base.string(S::T where T<:Field.ReceptiveField; sigdigits=2)
        M = ["$k=$(round(v;sigdigits))" for (k,v) in S.metrics]
        join(M, " ")
    end
    function Table.to_dataframe(F::ReceptiveField, pos...; kws...)
        F = Utils.to_dict(F)
        grid  = pop!(F, :grid_centers)
        props = pop!(F, :props)
        F = NamedTuple(F)
        if hasproperty(kws, :key_name)
            push!(kws.key_name, "field_prop")
        end
        Table.to_dataframe(F, pos...; grid=grid, props=props, kws...)
    end

    # Collection types
    abstract type FieldCollection end
    abstract type FieldDict  <: FieldCollection end
    abstract type FieldArray <: FieldCollection end
    sparse_to_full(G::Vector...) = [g for g in zip(ndgrid(G...)...)]
    sparse_to_full(G::Tuple{Vector}) = [g for g in zip(ndgrid(G...)...)]
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

    """
        resolution_to_width

    converts resolution to width, ie the number of points spread over a
    range of points into the inter-sample width
    """
    function resolution_to_width(resolution::OrderedDict,
            boundary::OrderedDict)::OrderedDict
        width = OrderedDict{String,Any}()
        for prop in keys(resolution)
            width[prop] = boundary[prop] / resolution[prop]
        end
        width
    end

    """
        get_boundary

    gets the boundary of each property
    """
    function get_boundary(behavior::DataFrame, props::Vector)::OrderedDict
        boundary   = OrderedDict{eltype(props),Any}()
        for prop in props
            boundary[prop]   = extrema(Utils.skipnan(behavior[!,prop]))
        end
        boundary
    end
    function center_to_edge(grid::AbstractVector)
        grid = collect(grid)
        Δ = median(diff(grid))
        δ = Δ/2
        grid = collect(minimum(grid)-δ:Δ:maximum(grid)+δ)
    end
    function edge_to_center(grid::AbstractArray)
        grid = collect(grid)
        grid = dropdims(mean([vec(grid[1:end-1]) vec(grid[2:end])], dims=2), dims=2)
    end

    """
        return_vals

    return a value matrix/dataframe for the requested
    properties in the DataFrame X
    """
    function return_vals(X::DataFrame, props::Vector)::Union{Vector{Float32},
                                                             Matrix{Float32}}
        Y = dropmissing(X[!, props])
        Float32.(hcat([Y[!,prop] for prop in props]...))
    end
    
    function isminutes(X::DataFrame)
        Utils.dextrema(X.time)[1] < 1440.0 # assumes less than 24 hour recording
    end

    function ensureTimescale!(X::DataFrame; kws...)
        if isminutes(X; kws...)
            transform!(X, :time => (x->x.*60) => :time)
        else
            X
        end
    end


    # Field-related submodules
    using Reexport
    include(srcdir("Field","operation.jl"))
    @reexport using .operation
    include(srcdir("Field","model.jl"))
    @reexport using .model
    include(srcdir("Field","fit.jl"))
    @reexport using .fit
    include(srcdir("Field","metrics.jl"))
    import .metrics
    include(srcdir("Field","recon.jl"))
    import .recon
    include(srcdir("Field","recon_process.jl"))
    import .recon_process
    include(srcdir("Field","adaptive.jl"))
    @reexport using .adaptive
    include(srcdir("Field","fixed.jl"))
    @reexport using .fixed
    include(srcdir("Field","preset.jl"))
    @reexport using .preset
    include(srcdir("Field","coactivity.jl"))
    @reexport using .coactivity

    """
        Base.:∈

    Method for asking if a point is inside a place field
    """
    function Base.:∈(point::Vector, field::ReceptiveField; hull=1)
        element(Singleton(point)) ∈ field.metrics[:hullseg_grid][hull]
    end

    #include(srcdir("Field","legacy.jl"))
    #import .legacy
end
