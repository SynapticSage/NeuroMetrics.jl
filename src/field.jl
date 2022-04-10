module field

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

    # Goal Vector Libraries
    #using DrWatson
    include("raw.jl")
    include("utils/SearchSortedNearest.jl/src/SearchSortedNearest.jl")
    include("utils.jl")
    # Submodules
    include("field/model.jl")
    export model
    include("field/info.jl")
    export info
    include("field/operation.jl")
    export operation

    rateConversion = 30
    export rateConversion

    function group(fields::Union{Tuple,Vector},
            names::Union{Vector{String},Vector{Symbol}};
            grouptype::Type=Dict, requireAll=false)

        keyset = keys(fields[1])
        for i in range(2,length(fields))
            keyset = union(keyset, keys(fields[i]))
        end
        keyset = Tuple(keyset)
        Grouping = Dict{typeof(keyset[1]), Any}()
        for key in keyset
            if requireAll
                notin = [!(key in keys(field)) for (name,field) in
                         zip(names,fields)]
                if any(notin)
                    continue
                end
            end
            if grouptype == Dict
                Grouping[key] = Dict(name=>field[key]
                                      for (name,field) in zip(names, fields)
                                      if key in keys(field))
            elseif grouptype == NamedTuple
                names = Symbol.(names)
                Grouping[key] = (;((name,field[key]) for (name,field) in
                                   zip(names, fields) if key in keys(field))...)
            else
                Grouping[key] = (grouptype)(field[key]
                                      for (name,field) in zip(names, fields)
                                      if key in keys(field))
            end
        end
        return Grouping
    end

    function getSettings(thing, props;
            resolution::Union{Vector{Int},Int,Nothing}=100,
            settingType::String="hist")

        if resolution isa Int
            resolution = repeat([resolution]; inner=length(props))
        end

        grid = OrderedDict{String,Any}()
        thing = dropmissing(thing);
        for prop in props
            grid[prop] = extrema(utils.skipnan(thing[!, prop]));
        end
        range_func_hist(start, stop, i) = collect(start : (stop-start)/resolution[i] : stop);
        range_func_kde(x,y,i) = x:((y-x)/resolution[i]):y
        edge_end(x,y,i) = range_func_kde(x,y,i)[begin:end-1]
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
            filters::Union{Dict,Nothing}=nothing, 
            props::Vector{String}=[],
            behfilter::Union{Dict, Bool}=nothing)

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
            println("augmented_props=$augmented_props")
            tmp, data = raw.filterAndRegister(beh, data; filters=filters,
                                              transfer=augmented_props,
                                              on="time")
            @assert length(unique(data.velVec)) > 1 "Fuck 2"
            if behfilter isa Bool
                if behfilter
                    println("Replacing behavior with filtered behavior")
                    beh = tmp;
                else
                    println("NOT replacing behavior with filtered behavior")
                end
            elseif behfilter isa Dict
                    println("replacing behavior with custom filter")
                beh = raw.filterTables(beh; filters=behfilter,
                                       lookupcols=nothing)[1]
            end
        end
        #@assert length(unique(data.velVec)) > 1 "Fuck 3"
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
            behfilter::Union{Bool, Dict}=true,
            filters::Union{Dict,Nothing}=nothing)

        # Add revelent columns to data and filter?
        beh, data = _handle_missing_cols_and_filters(copy(beh), copy(data);
                                                     filters=filters,
                                                     behfilter=behfilter,
                                                     props=props)

        if isempty(beh) || isempty(data)
            println("size(beh)=$(size(beh))")
            println("size(data)=$(size(data))")
            throw(ArgumentError("Your filtration leads to empty data"))
        end

        # HISTOGRAM
        H = (;hist=nothing, grid=nothing, occ=nothing, occzeroinds=nothing)
        if dohist
            println("props=$props")
            H = hist.fields(data, beh; props=props,
                                        savemem=savemem,
                                        resolution=resolution, 
                                        splitby=splitby,
                                        gaussian=gaussian);
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

    #Field = NamedTuple{(:hist, :kde, :cgrid, :egrid, :gridh, :gridk, :beh, :behdens), 
    #                   Tuple{Dict, Dict, Tuple, Tuple, Tuple, Vector, Array, Array}}

    function center_to_edge(grid)
        grid = collect(grid)
        Δ = median(diff(grid))
        δ = Δ/2
        grid = minimum(grid)-δ:Δ:maximum(grid)+δ
    end
    function edge_to_center(grid)
        grid = collect(grid)
        grid = dropdims(mean([vec(grid[1:end-1]) vec(grid[2:end])], dims=2), dims=2)
    end
    function to_density(field)
        field = field./nansum(vec(field))
    end
    function to_density(field, density)
        throw(InvalidStateException("to density not defined for this case yet"))
    end

    module hist
        using ..field
        import ..utils
        using DataFrames
        using StatsBase

        function h2d(thing::DataFrame, props::Vector{String}; grid=(),
                hist2dkws=Dict())
            thing = dropmissing(thing);
            P = tuple([convert(Vector{Float64}, x) for x in eachcol(thing[:, props])]...)
            return fit(Histogram, P, grid; hist2dkws...)
        end

        function fields(data::DataFrame, beh::DataFrame;
                splitby::Union{Nothing,Vector{Symbol},Vector{String},String}="unit",
                resolution::Union{Vector{Int},Int}=50, 
                savemem::Bool=true,
                props::Vector{String}=["x","y"],
                gaussian::Real = 0.0
            )
            if gaussian isa Int
                gaussian = convert(Float64, gaussian);
            end

            grid = field.getSettings(beh, props, 
                                    resolution=resolution,
                                    settingType="hist");
            println("Grid=$grid")
            println("resolution=$resolution")

            behDist = field.hist.h2d(beh, props, grid=grid);
            
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
                for i ∈ 1:length(groups)
                    if savemem
                        spikeDist[ugroups[i]] = Float32.(field_of_group(i).weights)
                    else
                        spikeDist[ugroups[i]] = field_of_group(i).weights
                    end
                end
            else
                if splitby != nothing
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
    export hist

    module kerneldens
        using ..field
        using DataFrames
        using StatsBase, KernelDensity
        using KernelDensitySJ
        #using DrWatson
        include("utils.jl")
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
                if bandwidths == nothing
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
                area_hist = sum(utils.skipnan(hist))
                area_kde  = sum(utils.skipnan(kde))
                kde = (area_hist/area_kde) .* kde;
            elseif kde isa Dict
                for key in keys(kde)
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
            grid = getSettings(beh, props, 
                                     resolution=resolution, settingType="kde")
            grid_hist = getSettings(beh, props, resolution=resolution,
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
                spikeDist = Dict(ugroups[i] => field_of_group(i) 
                                 for i ∈ 1:length(groups))
            else
                if splitby != nothing
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
                for i ∈ keys(spikeDist)
                    if spikeDist[i] == nothing
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
                @assert dist != nothing
                dist = pdf(spikeDist, grid...)
            end

            return (kde=dist, grid=grid, occ=behDist, occzeroinds=behzeroinds)
        end
    end
    export kerneldens

    """
    `plot`

    # proposed structure

    `show_fieldgroups` : plot a group of field collections, who are indexed by group keys

    `show_fields` : plot a collection of fields indexed by keys

    `show_field`  : plot a single field and its key

    """
    module plot
        using ..field
        using Plots, LaTeXStrings, Measures
        using Statistics
        using ProgressMeter
        import ..utils

        function transform_key(key; splitter::String="\n")
            tk(key) = replace(string(key), "("=>"",",)"=>"",", "=>splitter,
                                         ")"=>"", "\""=>"", " "=>"")
            if key isa Nothing
                annotation = ""
            elseif key isa String
                annotation = key
            else
                annotation = tk(key)
            end
        end

        function show_fieldgroups(group; groupgrid="row", as=Plots.plot,
                show_field_kws...)
            if groupgrid == nothing
                groupgrid = "rowcol"
            end
            if groupgrid isa String
                key = Tuple(keys(group))[1]
                N = length(group)
                if groupgrid == "row"
                    groupgrid = (N, 1)
                elseif groupgrid == "col"
                    groupgrid = (1, N)
                elseif groupgrid == "rowcol"
                    groupgrid = (ceil(sqrt(N)), floor(sqrt(N)))
                end
            end
            keys_ = Tuple(keys(group))
            N = length(keys_)
            if as == Plots.plot
                extrakws = (;nplots_factor=N, show_field_kws...)
            else
                extrakws = (; show_field_kws...,)
            end
            plots = [show_fields(group[key]; keyappend=key, extrakws...)
                     for key in keys_]
            if as == Plots.plot
                grid = [p.layout for p in plots]
                print(size(grid))
                grid = reshape(grid, groupgrid)
                plots=Plots.plot(plots..., grid=grid)
            elseif as == Dict || as == NamedTuple
                plots = as(zip(keys_,plots))
            else
                plots = as(plots)
            end
            return plots
        end
        function save_dict_fieldplots(plots, name; ext="pdf")
            @showprogress for (key, plot) in plots
                keyname = name * "_" * transform_key(key, splitter=",")
                if ext isa String
                    ext = [ext]
                end
                for e in ext
                    savefig(plot, keyname * ".$e")
                end
            end
        end
        function selectbackground()
            ndim = ndims(operation.selectind(F, 1))
        end
        function show_fields(F::Dict; fontscale=true, background=:grey30,
                textcolor=:white, as::Union{Type,<:Function}=Plots.plot,
                plotkws::NamedTuple=NamedTuple(), kws...)
            plotkws = (;selectbackground(F)..., plotkws...)
            if fontscale
                kws2 = (;nplots=length(F))
            else
                kws2 = ()
            end
            if as == Dict
                obj = as(key=>show_field(f; kws2..., key=key, kws...) for (key,f) in F)
            elseif as == Plots.plot
                plotkws=(;margin=-2mm, background_color=background,
                         foreground_color=background,
                         background_color_outside=background, plotkws...)
                kws = (;textcolor=textcolor, kws...)
                obj =  as((show_field(f; kws2..., key=key, kws...) for (key,f)
                           in sort(collect(F), by=x->x[1]))...; plotkws...)
            else
                obj =  as(show_field(f; kws2..., key=key, kws...) for (key,f) in F)
            end
            return obj
        end
        function show_field(FF::AbstractVector; func=Plots.bar,
                key::Union{String,NamedTuple,Nothing}=NamedTuple(),
                keyappend::Union{String,NamedTuple,Nothing}=nothing,
                keyprepend::Union{String,NamedTuple,Nothing}=nothing,
                grid=[], textcolor=:black, justification::Symbol=:bottom,
                fontsize=12, location=(1, 0.03), nplots=1, nplots_factor=1,
                quant::Vector{Float64}=[0.05,0.99], kws...)
            if all(isnan.(FF)) || isempty(FF)
                return Plots.plot()
            end
            ylims = quantile(utils.skipnan(vec(FF)), quant)
            if grid == []; extrakwargs = (xticks=[], yticks=[])
            else; extrakwargs = ()
            end
            kws = (;extrakwargs..., kws...)
            l = func(grid..., FF; kws...)
            annotate_field(l; key=key, keyappend=keyappend, keyprepend=keyprepend,
                           grid=grid, textcolor=textcolor,
                           justification=justification,
                           location=location,nplots=nplots,
                           nplots_factor=nplots_factor)
        end
        function show_field(FF::AbstractMatrix; 
                key::Union{String,NamedTuple,Nothing}=NamedTuple(),
                keyappend::Union{String,NamedTuple,Nothing}=nothing,
                keyprepend::Union{String,NamedTuple,Nothing}=nothing,
                grid=[], textcolor=:black, justification::Symbol=:bottom,
                fontsize=12, location=(1, 0.03), nplots=1, nplots_factor=1,
                quant::Vector{Float64}=[0.05,0.99], kws...)
            if all(isnan.(FF)) || isempty(FF)
                return Plots.plot()
            end
            clims = quantile(utils.skipnan(vec(FF)), quant)
            if grid == []; 
                extrakwargs = (xticks=[], yticks=[])
                grid = (); pos = ();
            else; extrakwargs = (); pos = grid;
            end
            kws = (;extrakwargs..., kws...)
            # ------
            # HEATMAP
            # ------
            hm = heatmap(pos..., FF; clims=Tuple(clims), 
                         colorbar=false, padding=(0,0),
                         c=cgrad(:acton,rev=false), showaxis=:no, kws...)
            #old cms = :linear_kryw_5_100_c67_n256
            #:lajolla
            annotate_field(hm; key=key, keyappend=keyappend, keyprepend=keyprepend,
                           grid=grid, textcolor=textcolor,
                           justification=justification,
                           location=location,nplots=nplots,
                           nplots_factor=nplots_factor)
        end
        function show_field(FF::AbstractArray{T,3}; 
                    key::Union{String,NamedTuple,Nothing}=NamedTuple(),
                    keyappend::Union{String,NamedTuple,Nothing}=nothing,
                    keyprepend::Union{String,NamedTuple,Nothing}=nothing,
                    grid=[], textcolor=:black, justification::Symbol=:bottom,
                    fontsize=12, location=(1, 0.03), nplots=1, nplots_factor=1,
                    quant::Vector{Float64}=[0.05,0.99], seriestype=:volume,
                    slice_dim=3, kws...) where T <: Real
            if all(isnan.(FF)) || isempty(FF)
                return Plots.plot()
            end
            if seriestype isa Symbol
                println("Volume path")
                kws = (; margins=-2mm, padding=(-1,-1), kws...)
                    #FF = replace(FF, NaN=>0)
                    t=text_annotate(key=key, keyappend=keyappend, 
                                    keyprepend=keyprepend,
                                   grid=grid, textcolor=textcolor,
                                   justification=justification,
                                   location=location, nplots=nplots,
                                   nplots_factor=nplots_factor)
                    kws = (; t..., kws...)
                    p = Plots.plot3d(FF; seriestype=seriestype, kws...)
                try
                    println("volume plotted")
                    println("annotation plotted")
                catch
                    @warn "Failed on key=$key"
                end
            else
                p = mapslices(x->show_field(x;key=key, keyappend=keyappend,
                                        keyprepend=keyprepend, grid=grid,
                                        textcolor=textcolor,
                                        justification=justification,
                                        fontsize=fontsize, location=location,
                                        nplots=nplots,
                                        nplots_factor=nplots_factor,
                                        quant=quant),
                          FF, dims=slice_dim)
            end
            return p
        end
        function annotate_field(hm;
                    key::Union{String,NamedTuple,Nothing}=NamedTuple(),
                    keyappend::Union{String,NamedTuple,Nothing}=nothing,
                    keyprepend::Union{String,NamedTuple,Nothing}=nothing,
                    grid=[], textcolor=:black, justification::Symbol=:bottom,
                    fontsize=12, location=(1, 0.03), nplots=1, nplots_factor=1)
            # ANNOTATION OF KEY
            # -----------------
            fontsize = max(Int(round(fontsize/(nplots*nplots_factor))),3)
            #println("fontsize=$fontsize")
            get_loc(p,l) = p[1] + (p[2]-p[1])*l
            x = get_loc(xlims(), location[1])
            y = get_loc(ylims(), location[2])
            color = textcolor
            font = :bold
            annotation = [transform_key(key)]
            if keyappend != nothing
                push!(annotation, transform_key(keyappend))
            end
            if keyprepend != nothing
                pushfirst!(annotation, transform_key(keyappend))
            end
            annotation = join(annotation, "\n")
            t = text(annotation, justification, color, :right, :bold, pointsize=fontsize)
            annotate!(hm, x, y, t) 
        end
        function text_annotate(;
                    key::Union{String,NamedTuple,Nothing}=NamedTuple(),
                    keyappend::Union{String,NamedTuple,Nothing}=nothing,
                    keyprepend::Union{String,NamedTuple,Nothing}=nothing,
                    grid=[], textcolor=:black, justification::Symbol=:bottom,
                    fontsize=12, location=(1, 0.03), nplots=1, nplots_factor=1)
            # ANNOTATION OF KEY
            # -----------------
            fontsize = max(Int(round(fontsize/(nplots*nplots_factor))),3)
            #println("fontsize=$fontsize")
            get_loc(p,l) = p[1] + (p[2]-p[1])*l
            x = get_loc(xlims(), location[1])
            y = get_loc(ylims(), location[2])
            color = textcolor
            font = :bold
            annotation = [transform_key(key)]
            if keyappend != nothing
                push!(annotation, transform_key(keyappend))
            end
            if keyprepend != nothing
                pushfirst!(annotation, transform_key(keyappend))
            end
            annotation = join(annotation, "\n")
            return (;title=annotation, titlefontsize=fontsize,
                    titlefontcolor=textcolor, titlefontvalign=justification)
        end
        function title_annotate(hm; kws...)
            annotation = text_annotate(;kws...).title
            hm.attr[:plot_title]    = annotation
            hm.attr[:plot_titlefontvaling] = justification
            hm.attr[:plot_titlefontcolor] = color
            hm.attr[:window_title] = annotation
        end
    end
    
    export plot

end

