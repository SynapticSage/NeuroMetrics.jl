module legacy

    using DataFrames
    using DataStructures: OrderedDict
    import GeometricalPredicates

    import DIutils
    import Load

    # ------------- LEGACY CODE ---------------------------------------------
    function getSettings(thing, props;
            resolution::Union{Vector{Int},Int,Nothing}=100,
            settingType::String="hist")

        if resolution isa Int
            resolution = repeat([resolution]; inner=length(props))
        end

        grid = OrderedDict{String,Any}()
        thing = dropmissing(thing);
        for prop in props
            grid[prop] = extrema(DIutils.skipnan(thing[!, prop]));
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

        user_specified_filter = (filters !== nothing)
        missing_columns_in_df = !all(in.(props,[names(data)]))

        if user_specified_filter || missing_columns_in_df
            if user_specified_filter
                filter_props = [key for (key,_) in filters
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

    ### Input
    `beh`::DataFrame
    `data`::DataFrame
    `resolution`     -- ::Int=40 Resolution to sample our fields at
    `hist2kde_ratio` -- ::Real=40/100Ratio of hist to kde sampling resolution
    `props`          -- ::Vector{String}=["x","y"]Which props to make a field of
    `gaussian`       -- ::Real=0 Gaussian kernel? If 0, none applied (kde doesn't need it)
    `dokde`          -- ::Bool=true Do KDE at all?
    `dohist`         -- ::Bool=true Do hist at all?
    `normkde`        -- ::Bool=true Normalize kde area by hist area (make them occupy same area)
    `behfilter`      -- ::Bool=true If false, we skip applying the filter to behavior (can be useful for normalizing purposes, because non-occupied areas are naned)
    `filters`        -- ::Union{Dict,Nothing}=nothing
    ### Output
    `X` -- ::NamedTuple
    """
    function get_fields(beh::DataFrame, data::DataFrame; 
            resolution::Union{Int, Vector{Int}, Tuple{Int}}=40, 
            hist2kde_ratio::Real=1,
            props::Vector{String}=["x","y"],
            splitby::Union{Vector{Symbol},Vector{String},Nothing,String}="unit",
            gaussian::Real=0,
            dokde::Bool=false,
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
        if boundary !== nothing
            H  =_enforce_boundary(H, boundary)
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

    function to_density(field::AbstractArray)
        field = field./nansum(vec(field))
    end
    function to_density(field, density)
        @error "not defined yet"
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
        points = __get_point.(Iterators.product(grid ...))
        boundary = [__get_point(r) for r in eachrow(boundary)]
        polygon = GeometricalPredicates.Polygon2D(boundary...)
        answers = (!).(GeometricalPredicates.inpolygon.([polygon], points))
        H[answers] .= NaN
    end

    module hist

        using ..legacy
        import DIutils

        using DataFrames
        using StatsBase: fit
        using StatsBase
        using Infiltrator

        function h2d(thing::DataFrame, props::Vector{String}; grid=(),
                hist2dkws=Dict())
            thing = dropmissing(thing);
            P = tuple([convert(Vector{Float64}, x) for x in
                       eachcol(thing[:, props])]...)
            return fit(Histogram, P, grid; hist2dkws...)
        end

        function fields(data::DataFrame, beh::DataFrame;
                splitby::Union{Nothing,Vector{Symbol},Vector{String},String}="unit",
                resolution::Union{Vector{Int},Int}=50, 
                savemem::Bool=true,
                props::Vector{String}=["x","y"],
                gaussian::Real=0.0
            )
            if gaussian isa Int
                gaussian = convert(Float64, gaussian);
            end

            grid = legacy.getSettings(beh, props, 
                                    resolution=resolution,
                                    settingType="hist");
            @debug "Grid=$grid"
            @debug "resolution=$resolution"

            behDist = hist.h2d(beh, props, grid=grid);
            
            ith_group(i) = DataFrame(groups[i])
            field_of_group(i) = hist.h2d(ith_group(i), props, grid=grid)
            if splitby isa Vector{Symbol}
                splitby = String.(splitby)
            end
            if all(in.(splitby, [names(data)]))
                groups = groupby(data, splitby)
                ugroups =
                collect(zip(sort(collect(groups.keymap),by=x->x[2])...))[1]
                ugroups = 
                [(;zip(groups.cols,ugroup)...) for ugroup in ugroups]
                #spikeDist = Dict(ugroups[i] => field_of_group(i) 
                #                 for i ∈ 1:length(groups))
                spikeDist = Dict{typeof(ugroups[1]),Any}();
                @inbounds for i ∈ 1:length(groups)
                    if savemem
                        spikeDist[ugroups[i]] = Float32.(field_of_group(i).weights)
                    else
                        spikeDist[ugroups[i]] = field_of_group(i).weights
                    end
                end
            else
                if splitby !== nothing
                    @warn "splitby not empty, yet not all columns in df"
                end
                spikeDist = hist.h2d(data, props, grid=grid)
            end

            # Transform behavior
            behDist_nanWhere0 = copy(behDist.weights);
            if savemem
                behDist_nanWhere0 = convert.(Float32, behDist_nanWhere0);
            else
                behDist_nanWhere0 = convert.(Float64, behDist_nanWhere0);
            end
            behzeroinds = behDist_nanWhere0 .== 0
            behDist_nanWhere0[behzeroinds] .= NaN;
            behDist.weights[behDist.weights .== 0] .= 1;
            
            return (hist=spikeDist, grid=grid, occ=behDist_nanWhere0,
                    occzeroinds=behzeroinds)
        end
    end

    module kerneldens

        using ..legacy
        using DataFrames
        using StatsBase, KernelDensity
        using KernelDensitySJ
        #using Infiltrator
        #using DrWatson
        import DIutils

        function KDE(data, props; bandwidth=:silverman)
            data = dropmissing(data[:,props])
            #for col in props
            #    inds = data[:, col] .== NaN
            #    data[inds, col] .= 0
            #end
            if isempty(data) || data isa Nothing
                return nothing
            else
                fail=false
                bandwidths=nothing
                if bandwidth == :sj
                    try
                        bandwidths = tuple([bwsj(data[:,prop]) for prop in props]...)
                    catch e
                        fail=true
                    end
                    if fail || !all(bandwidths .> 0) 
                        return nothing
                    end
                elseif bandwidth == :silverman
                    bandwidths = nothing
                else
                    bandwidths = bandwidth
                end
                if bandwidths === nothing
                    return kde(tuple([convert(Vector{Float32}, x) 
                                      for x in eachcol(data[:, props])]...
                                    )
                              )
                else
                    return kde(tuple([convert(Vector{Float32}, x) 
                                      for x in eachcol(data[:, props])]...
                                    ); bandwidth=bandwidths
                              )
                end
            end
        end


        function gridpdf(kde, grid)
            grid = [grid[prop] for prop in keys(grid)]
            pdf(kde, grid...)
        end

        """
            norm_kde_by_histcount

        normalizes kernel density estimates to have the same number of event
        counts as the histogram, so that it will accurately represent peak
        firing/ripple rates
        """
        function norm_kde_by_histcount(kde::Union{AbstractArray, Dict}, 
                hist::Union{AbstractArray,Dict}, skip_keyerror=false)
            #if typeof(hist) ≠ typeof(kde)
            #    throw(ArgumentError("type of hist should match type of kde\n...type(hist)=$(supertype(typeof(hist))) != $(supertype(typeof(kde)))"))
            #end
            if kde isa AbstractArray
                area_hist = sum(DIutils.skipnan(hist))
                area_kde  = sum(DIutils.skipnan(kde))
                kde = (area_hist/area_kde) .* kde;
            elseif kde isa Dict
                @inbounds for key in keys(kde)
                    try
                        kde[key] = norm_kde_by_histcount(kde[key], hist[key])
                    catch KeyError
                        if !skip_keyerror
                            throw(KeyError("Key=$key missing from hist"))
                        end
                    end
                end
            else
                throw(ArgumentError)
            end
            return kde
        end

        function fields(data::DataFrame, beh::DataFrame;
                splitby::Union{Nothing,Vector{String},Vector{Symbol},String}="unit",
                resolution::Union{Int,Vector{Int}}=200,
                savemem::Bool=true,
                props::Vector{String}=["x","y"])

            # Grid settings
            grid = legacy.getSettings(beh, props, 
                                     resolution=resolution, settingType="kde")
            #@infiltrate
            grid_hist = legacy.getSettings(beh, props, resolution=resolution,
                                    settingType="hist")
            # Behavioral distribution
            behDist = KDE(beh, props);
            behDist_hist = hist.h2d(beh, props, grid=grid_hist).weights;

            # Data distribution applying splits if applicable
            ith_group(i) = DataFrame(groups[i])
            field_of_group(i) = KDE(ith_group(i), props)
            if splitby isa Vector{Symbol}
                splitby = String.(splitby)
            end

            if all(in.(splitby, [names(data)]))
                groups = groupby(data, splitby)
                ugroups = collect(zip(sort(collect(groups.keymap),
                                           by=x->x[2])...))[1]
                ugroups = 
                [(;zip(groups.cols,ugroup)...) for ugroup in ugroups]
                spikeDist = Dict(ugroups[i] => try field_of_group(i) catch nothing end
                                 for i ∈ 1:length(groups)) # TODO replace nothing with array of nan
            else
                if splitby !== nothing
                    @warn "splitby not empty, yet not all columns in df"
                end
                spikeDist = KDE(data[!,:], props)
            end

            # Occupancy
            behDist = pdf(behDist, grid...);
            zero_fraction = 0.0000000001
            behDist[behDist .<= zero_fraction] .= 1
            behzeroinds = behDist_hist .<= 0
            if savemem
                behDist = Float32.(behDist)
            end

            # Spike count
            if spikeDist isa Dict
                dist = Dict{typeof(ugroups[1]),Any}();
                @inbounds for i ∈ keys(spikeDist)
                    if spikeDist[i] === nothing
                        dist[i] = nothing
                    else
                        if savemem
                            dist[i] = Float32.(pdf(spikeDist[i], grid...))
                        else
                            dist[i] = pdf(spikeDist[i], grid...)
                        end
                    end
                end
            else
                @assert dist !== nothing
                dist = pdf(spikeDist, grid...)
            end

            return (kde=dist, grid=grid, occ=behDist, occzeroinds=behzeroinds)
        end

    end

    

end
