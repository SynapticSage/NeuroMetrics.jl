module Field

    export getSettings
    export get_fields, Field
    export skipnan

    # Julia Packages
    using DataFrames
    using ImageFiltering
    import DataStructures: OrderedDict
    using DataStructures
    using Statistics
    using NaNStatistics
    using Infiltrator
    using GeometricalPredicates: inpolygon

    # Goal Vector Libraries
    using DrWatson
    import Load
    import Utils
    

    rateConversion = 30
    export rateConversion

    function getSettings(thing, props;
            resolution::Union{Vector{Int},Int,Nothing}=100,
            settingType::String="hist")

        if resolution isa Int
            resolution = repeat([resolution]; inner=length(props))
        end

        grid = OrderedDict{String,Any}()
        thing = dropmissing(thing);
        for prop in props
            grid[prop] = extrema(Utils.skipnan(thing[!, prop]));
        end

        range_func_hist(start, stop, i) = collect(start : (stop-start)/resolution[i] : stop)
        range_func_kde(x,y,i) = x:((y-x)/resolution[i]):y
        edge_end(x,y,i)   = range_func_kde(x,y,i)[begin:end-1]
        edge_start(x,y,i) = range_func_kde(x,y,i)[begin+1:end]
        mp(x,y,i) = dropdims(diff(hcat(edge_start(x,y,i),edge_end(x,y,i)),
                                   dims=2); dims=2) +  edge_start(x,y,i);

        get_extrema(prop) = grid[prop][1], grid[prop][2]

        if lowercase(settingType) == "hist"
            grid    = OrderedDict(x[1] => range_func_hist(x[2]..., i) for (i,x) in
                           Iterators.enumerate(grid))
            grid    = Tuple(grid[name] for name in props)
        elseif lowercase(settingType) == "kde"
            grid = [mp(get_extrema(prop)...,i ) for (i,prop) in Iterators.enumerate(props)]
        end

        return grid
    end

    function _handle_missing_cols_and_filters(beh::DataFrame, 
            data::DataFrame; 
            filters::Union{AbstractDict,Nothing}=nothing, 
            props::Vector{String}=[],
            behfilter::Union{AbstractDict, Bool}=nothing)

        user_specified_filter = (filters != nothing)
        missing_columns_in_df = !all(in.(props,[names(data)]))

        if user_specified_filter || missing_columns_in_df
            if user_specified_filter
                filter_props = [key for (key,value) in filters
                                if key != All()]
            else
                filter_props = ()
            end
            augmented_props = [filter_props..., props...]
            @debug "augmented_props=$augmented_props"
            tmp, data = Load.filterAndRegister(beh, data; filters=filters,
                                              transfer=augmented_props,
                                              on="time")
            @assert length(unique(data.velVec)) > 1 "Empty data!"
            if behfilter isa Bool
                if behfilter
                    @debug "Replacing behavior with filtered behavior"
                    beh = tmp;
                else
                    @debug "NOT replacing behavior with filtered behavior"
                end
            elseif behfilter isa Dict
                    @debug "replacing behavior with custom filter"
                beh = Load.filterTables(beh; filters=behfilter,
                                       lookupcols=nothing)[1]
            end
        end
        beh, data
    end

    """
        get_fields

    distills the process of getting fields into one interface

    # Arguments
        beh::DataFrame

        data::DataFrame

        resolution::Int=40 
    Resolution to sample our fields at

        hist2kde_ratio::Real=40/100,
    Ratio of hist to kde sampling resolution

        props::Vector{String}=["x","y"],
    Which props to make a field of

        gaussian::Real=0,
    Gaussian kernel? If 0, none applied (kde doesn't need it)

        dokde::Bool=true,
    Do KDE at all?

        dohist::Bool=true,
    Do hist at all?

        normkde::Bool=true,
    Normalize kde area by hist area (make them occupy same area)

        behfilter::Bool=true
    If false, we skip applying the filter to behavior (can be useful
    for normalizing purposes, because non-occupied areas are naned)

        filters::Union{Dict,Nothing}=nothing)
    # Output
        X::named_tuple
    """
    function get_fields(beh::DataFrame, data::DataFrame; 
            resolution::Union{Int, Vector{Int}, Tuple{Int}}=40, 
            hist2kde_ratio::Real=1,
            props::Vector{String}=["x","y"],
            splitby::Union{Vector{Symbol},Vector{String},Nothing,String}="unit",
            gaussian::Real=0,
            dokde::Bool=true,
            dohist::Bool=true,
            normkde::Bool=true,
            savemem::Bool=true,
            boundary::Union{Nothing, DataFrame}=nothing,
            behfilter::Union{Bool, AbstractDict}=true,
            filters::Union{AbstractDict,Nothing}=nothing)

        # Add revelent columns to data and filter?
        beh, data = _handle_missing_cols_and_filters(copy(beh), copy(data);
                                                     filters=filters,
                                                     behfilter=behfilter,
                                                     props=props)

        if isempty(beh) || isempty(data)
            @debug "size(beh)=$(size(beh))"
            @debug "size(data)=$(size(data))"
            throw(ArgumentError("Your filtration leads to empty data"))
        end

        # HISTOGRAM
        H = (;hist=nothing, grid=nothing, occ=nothing, occzeroinds=nothing)
        if dohist
            @debug "props=$props"
            H = hist.fields(data, beh; props=props,
                                        savemem=savemem,
                                        resolution=resolution, 
                                        splitby=splitby,
                                        gaussian=gaussian);
        end

        # Boundary?
        if boundary != nothing
            H  =_enforce_boundary(H)
        end

        # KERNEL DENSITY
        atmost2d = (length(props) <= 2)
        K = (;kde=nothing, grid=nothing, occ=nothing, occzeroinds=nothing)
        if dokde && atmost2d
            K = kerneldens.fields(data, beh; props=props,
                                         savemem=savemem,
                                         splitby=splitby,
                                         resolution=Int.(resolution ./
                                                        hist2kde_ratio));
            if normkde
                kdfields = kerneldens.norm_kde_by_histcount(K.kde, H.hist)
                K = (;kde=kdfields, grid=K.grid, occ=K.occ, occzeroinds=K.occzeroinds)
            end
        end

        gride = H.grid
        gridc = edge_to_center.(gride)
        
        out =  (Cₕ=H.hist, Cₖ=K.kde, occ=H.occ, occR=to_density(H.occ),
                occzeroinds=H.occzeroinds, cgrid=gridc, egrid=gride,
                gridh=H.grid, gridk=K.grid, dims=props)
        out = operation.occnorm(out)
        return out
    end

    function center_to_edge(grid::AbstractVector)
        grid = collect(grid)
        Δ = median(diff(grid))
        δ = Δ/2
        grid = minimum(grid)-δ:Δ:maximum(grid)+δ
    end
    function edge_to_center(grid::AbstractArray)
        grid = collect(grid)
        grid = dropdims(mean([vec(grid[1:end-1]) vec(grid[2:end])], dims=2), dims=2)
    end
    function to_density(field::AbstractArray)
        field = field./nansum(vec(field))
    end
    function to_density(field, density)
        throw(InvalidStateException("to density not defined for this case yet"))
    end

    function _enforce_boundary(H::NamedTuple, boundary::Matrix)
        @infiltrate
        for field in RF._dictfields
            h =getproperty(H, field)
            grid = H.cgrid
            h = _enforce_boundary(h, grid, boundary)
        end
    end
    function _enforce_boundary(H::Dict, grid::Tuple, boundary::Matrix)
        Dict(k => _enforce_boundary(v, grid, boundary) for (k,v) in H)
    end

    __get_point(x::Tuple) = GeometricalPredicates.Point(x...)
    __get_point(x::AbstractVector) = GeometricalPredicates.Point(x...)
    function _enforce_boundary(H::Array, grid::Tuple, boundary::Matrix)
        @infiltrate
        points = __get_point.(Iterators.product(grid ...))
        boundary = [__get_point(r) for r in eachrow(boundary)]
        polygon = GeometricalPredicates.Polygon2D(boundary...)
        answers = (!).(GeometricalPredicates.inpolygon.([polygon], points))
        H[answers] .= NaN
    end

    function get_RF()
    end

    # Field-related submodules
    using Reexport
    include(srcdir("Field","operation.jl"))
    @reexport using .operation
    include(srcdir("Field","model.jl"))
    @reexport using .model
    include(srcdir("Field","fit.jl"))
    @reexport using .fit
    include(srcdir("Field","plot.jl"))
    @reexport using .plot
    #include(srcdir("Field","utils.jl"))
    #@reexport using .utils
    include(srcdir("Field","hist.jl"))
    import .hist
    include(srcdir("Field","kerneldens.jl"))
    import .kerneldens
    include(srcdir("Field","info.jl"))
    import .info
    include(srcdir("Field","recon.jl"))
    import .recon
    include(srcdir("Field","recon_process.jl"))
    import .recon_process
    include(srcdir("Field","RF.jl"))
    import .RF
    include(srcdir("Field","adaptive.jl"))
    import .adaptive

end

