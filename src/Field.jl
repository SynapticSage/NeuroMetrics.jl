__precompile__(false)
module Field

    export get_fields
    export ReceptiveField

    # Julia Packages
    using DataFrames, ImageFiltering, DataStructures, Statistics, NaNStatistics, 
          Infiltrator, LazySets, Statistics
    import DataStructures: OrderedDict
    using GeometricalPredicates: inpolygon

    # Goal Vector Libraries
    using DrWatson
    using DIutils
    import DI
    import DIutils.binning: return_vals

    rateConversion = 30
    export rateConversion

	WIDTHS = OrderedDict(
	    "x"=>5f0, "y"=>5f0, "currentPathLength"=>4f0,
        "currentAngle"=>Float32(2pi/40),
        "stopWell"=>0.5f0
	)
    getwidths(props::Vector{String};scale=1)::Vector{Float32} = [Field.WIDTHS[prop]*scale for prop in props]
    getwidthsdict(props::Vector{String};scale=1)= Dict(prop=>Field.WIDTHS[prop]*scale for prop in props )


    # Our basic RECEPTIVE-FIELD type components
    abstract type ReceptiveField end
    ReceptiveFields = AbstractArray{<:ReceptiveField}

    function Base.push!(R::ReceptiveField, measurement_name::Symbol, measurement_value)
        push!(R.metrics, measurement_name => measurement_value)
    end
    function Base.string(S::T where T<:Field.ReceptiveField; sigdigits=2)
        M = ["$k=$(round(v; sigdigits))" for (k,v) in S.metrics]
        join(M, " ")
    end
    function Table.to_dataframe(F::ReceptiveField, pos...; kws...)
        F = DIutils.to_dict(F)
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


    function center_to_edge(grid::AbstractVector)
        grid = collect(grid)
        Δ = median(diff(grid))
        δ = Δ/2
        grid = collect(minimum(grid)-δ:Δ:maximum(grid)+δ)
    end
    function edge_to_center(edges::AbstractArray)
        grid = collect(edges)
        dropdims(mean([vec(grid[1:end-1]) vec(grid[2:end])], dims=2), dims=2)
    end
    
    function isminutes(df::DataFrame)
        DIutils.dextrema(df.time)[1] < 1440.0 # assumes less than 24 hour recording
    end

    function ensureTimescale!(df::DataFrame; kws...)
        if isminutes(df; kws...)
            transform!(df, :time => (x->x.*60) => :time)
        else
            df
        end
    end


    # Field-related submodules
    include("Field/fit.jl")
    include("Field/metrics.jl")
    include("Field/recon_process.jl")
    include("Field/adaptive.jl")
    #include(srcdir("Field","fixed.jl"))
    #@reexport using .fixed
    include("Field/preset.jl")
    include("Field/coactivity.jl")

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
