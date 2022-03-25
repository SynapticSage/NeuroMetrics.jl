module field

    export getSettings
    export get_fields
    export _occNormField
    export skipnan
    export field_to_dataframe, fields_to_dataframe

    # Julia Packages
    using DataFrames
    using ImageFiltering
    using LazyGrids: ndgrid
    import DataStructures: OrderedDict

    # Goal Vector Libraries
    using DrWatson
    @assert isfile(srcdir("raw.jl"))
    skipnan(x) = Iterators.filter(!isnan, x)
    include(srcdir("raw.jl"))
    include(srcdir("utils", "SearchSortedNearest.jl", "src",
                       "SearchSortedNearest.jl"))

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
                notin = [!(key in keys(field)) for (name,field) in zip(names,fields)]
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
            grid[prop] = extrema(skipnan(thing[!, prop]));
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
            behfilter::Union{Bool, Dict}=true,
            filters::Union{Dict,Nothing}=nothing)

        # Add revelent columns to data and filter?
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
            lookupcols = Dict((source=1,target=2) => augmented_props)
            tmp, data = raw.filterTables(beh, data; filters=filters,
                                              lookupcols=lookupcols)
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


        if isempty(beh) || isempty(data)
            throw(ArgumentError("Your filtration leads to empty data"))
        end

        if dohist
            println("props=$props")
            fields, gridf = hist.fields(data, beh; props=props,
                                        resolution=resolution, 
                                        splitby=splitby,
                                        gaussian=gaussian);
        else
            fields, gridf = nothing, nothing
        end
        atmost2d = (length(props) <= 2)
        if dokde && atmost2d
            kdfields, gridk = kerneldens.fields(data, beh; props=props,
                                         splitby=splitby,
                                         resolution=Int.(resolution ./
                                                        hist2kde_ratio));
            if normkde
                kdfields = kerneldens.norm_kde_by_histcount(kdfields, fields)
            end
        else
            kdfields, gridk = nothing, nothing
        end
        
        return (hist=fields, kde=kdfields, grid=gridk, gridf=gridf)
    end

    function _get_edge(grid)
        egrid
    end


    """
    +purpose: converts a field dictionary (multiple units) into a field dataframe
    """
    function fields_to_dataframe(fields; other_labels=Dict(),
            key_name::String="unit")

        D = DataFrame()
        for (key, field) in fields
            other_labels[key_name] = key
            append!(D, fields_to_dataframe(fields[key]; 
                                           other_labels=other_labels))
        end
        return D
    end


    """
        field_to_dataframe(field::AbstractArray; other_labels=Dict())

    +purpose: converts a single field matrix into a field dataframe
    """
    function field_to_dataframe(field::AbstractArray; other_labels=Dict())
        D = ndgrid((1:size(field,i) for i in 1:ndims(field))...)
        D = Dict{String,Any}("dim_$d"=>vec(D[d])
                              for d in 1:ndims(field))
        for (label, value) in other_labels
            print(label, value)
            D[label] = value
        end
        D["field"] = vec(field);
        D = DataFrame(D)
    end

    """
        fields_to_dataframe(field::AbstractArray; other_labels=Dict())

    +purpose: converts a nested dict of fields into a dataframe
    """
    function fields_to_dataframe(fields::Dict; other_labels=Dict())
    end

    module hist
        using ..field
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
                spikeDist = Dict(ugroups[i] => field_of_group(i) 
                                 for i ∈ 1:length(groups))
            else
                if splitby != nothing
                    @warn "splitby not empty, yet not all columns in df"
                end
                spikeDist = hist.h2d(data, props, grid=grid)
            end

            # Transform behavior
            behDist_nanWhere0 = copy(behDist.weights);
            behDist_nanWhere0 = convert.(Float64, behDist_nanWhere0);
            behzeroinds = behDist_nanWhere0 .== 0
            behDist_nanWhere0[behzeroinds] .= NaN;
            behDist.weights[behDist.weights .== 0] .= 1;
            
            if spikeDist isa Dict
                dist = Dict{typeof(ugroups[1]), Any}()
                for i ∈ keys(spikeDist)
                    dist[i] = _occNormField(spikeDist[i].weights,
                                                  behDist.weights;
                                                  behzeroinds=behzeroinds,
                                                  gaussian=gaussian)
                end
            else
                dist = _occNormField(spikeDist.weights, behDist.weights;
                                           behzeroinds=behzeroinds,
                                           gaussian=gaussian)
            end

            return dist, grid
        end
    end
    export hist

    module kerneldens
        using ..field
        using DataFrames
        using StatsBase, KernelDensity
        using KernelDensitySJ
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
            if typeof(hist) ≠ typeof(kde)
                throw(ArgumentError("type of hist should match type of kde"))
            end
            if kde isa AbstractArray
                area_hist = sum(skipnan(hist))
                area_kde  = sum(skipnan(kde))
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
                props::Vector{String}=["x","y"])

            # Grid settings
            grid = getSettings(beh, props, 
                                     resolution=resolution, settingType="kde")
            grid_hist = getSettings(beh, props, 
                                          resolution=resolution, settingType="hist")
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

            behDist = pdf(behDist, grid...);
            zero_fraction = 0.0000000001
            behDist[behDist .<= zero_fraction] .= 1
            behzeroinds = behDist_hist .<= 0

            if spikeDist isa Dict
                dist = Dict{typeof(ugroups[1]),Any}();
                for i ∈ keys(spikeDist)
                    if spikeDist[i] == nothing
                        dist[i] = nothing
                    else
                        cell = pdf(spikeDist[i], grid...)
                        dist[i] = _occNormField(cell, behDist; 
                                                      gaussian=0.0,
                                                      behzeroinds=behzeroinds)
                    end
                end
            else
                dist = pdf(spikeDist, grid...)
                @assert dist != nothing
                println(size(dist), typeof(dist))
                dist = _occNormField(dist, behDist; gaussian=0.0, 
                                           behzeroinds=behzeroinds)
            end
            return dist, grid
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
        function show_fields(F::Dict; fontscale=true,
                background=:grey30, textcolor=:white,
                as::Union{Type,<:Function}=Plots.plot,
                plotkws::NamedTuple=NamedTuple(), kws...)
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
                obj =  as((show_field(f; kws2..., key=key, kws...) for (key,f) in sort(collect(F), by=x->x[1]))...;
                         plotkws...)
            else
                obj =  as(show_field(f; kws2..., key=key, kws...) for (key,f) in F)
            end
            return obj
        end
        function show_field(F; key::Union{String,NamedTuple,Nothing}=NamedTuple(),
                keyappend::Union{String,NamedTuple,Nothing}=nothing,
                keyprepend::Union{String,NamedTuple,Nothing}=nothing,
                xy=[], textcolor=:black, justification::Symbol=:bottom,
                fontsize=12, location=(1, 0.03), nplots=1, nplots_factor=1,
                quant::Vector{Float64}=[0.05,0.99], kws...)
            if ((typeof(F) <: AbstractArray) == false) || 
                all(isnan.(F))
                return Plots.plot()
            end
            clims = quantile(skipnan(vec(F)), quant)
            if xy == []
                kwargs = (xticks=[], yticks=[])
            else
                kwargs = ()
            end
            # HEATMAP
            # ------
            hm = heatmap(xy..., F; clims=Tuple(clims), kwargs...,
                         colorbar=false, padding=(0,0),
                         c=cgrad(:acton,rev=false), showaxis=:no, kws...)
            #old cms = :linear_kryw_5_100_c67_n256
            #:lajolla
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
    end
    export plot

    module operation
        using Statistics
        using Bootstrap
        """
        field_operation

        used to collapse a field onto 1D
        """
        function unary(field, operation::Function=mean;
            kws...)
            if typeof(field) <: AbstractArray
                field = operation(field; kws...)
            elseif field isa Dict
                field = Dict(name=>unary(f, operation; kws...) for (name,f) in
                             field)
            end
            return field
        end
        """
        field_operation

        used to do an operation between fields
        """
        function binary(fieldA, fieldB, operation::Function=(./); kws...)
            if typeof(fieldA) <: AbstractArray
                field = operation(fieldA, fieldB; kws...)
            elseif fieldA isa Dict
                field = Dict(namesA=>binary(fieldA[namesA], fieldB[namesA],
                                          operation; kws...) for namesA in
                             keys(fieldA))
            else
                throw(TypeError("fielld is wrong type $(typeof(fieldA))"))
            end
            return field
        end
    end

    function _occNormField(data::AbstractArray,
            behDist::Union{Array,Matrix};
            behzeroinds::Union{AbstractArray,Nothing}=nothing,
            gaussian::Real=0.0)
        X = convert(Array{Float64}, data);
        X = replace(X, NaN=>0)
        #zeroinds = X .== 0
        X = X ./ behDist;
        if gaussian != 0.0
            X = imfilter(X, Kernel.gaussian(gaussian))
        end
        if behzeroinds != nothing
            X[behzeroinds] .= NaN;
        end
        return X
    end

    module model
        using DrWatson
        using StatsBase
        using ProgressMeter
        using DataFrames
        include(srcdir("utils", "SearchSortedNearest.jl", "src",
                       "SearchSortedNearest.jl"))

        """
        data

        gets a form that represents the data

        get the D of spikes per behavioral Δt
        """
        function data(spikes, behavior; grid=nothing, props=nothing,
                toleranceₑ=0.1, as="group", fixbigfact=true)
            spikes = copy(spikes)

            # Map each spike to its behavioral bin
            nearest_points = SearchSortedNearest.searchsortednearest
            spikes[!,"nearest"] = nearest_points.([behavior.time], spikes.time)

            # Filter out spikes that do not match the tolerance
            Δϵ = abs.(behavior.time[spikes.nearest] .- spikes.time)
            spikes = spikes[Δϵ .< toleranceₑ, :]

            # Split by tetrode and behavioral bin ⟶   get D!
            D = combine(groupby(spikes, [:unit, :nearest], sort=true), nrow)
            D = unstack(D, :unit, :nrow)
            for (c, name) in zip(1:size(D, 2), names(D))
                col = D[!, c]
                D[ismissing.(col), c] .= 0
                if fixbigfact && name != "nearest"
                    D[col .> 20, c] .= 20 # maximumal non-big() factorial
                end
            end

            # Let's lookup the linear index into a field for every behavioral
            # time
            if grid != nothing && props != nothing
                # Lookup each of the properties that we need to look into the grid fro
                properties = behavior[D[!, "nearest"], props]
                # From there, we have to find where each point in those proeprties
                # lives inside of the grid
                @assert length(grid) == length(props) "UH oh, look at the length of your passed in objects"
                properties = Matrix(properties)
                behavior[!, "fieldIndex"] = zeros(Int64, size(behavior,1))
                @showprogress for (i,row) in enumerate(eachrow(properties))
                    W = fit(Histogram, Tuple(([r] for r in row)), grid).weights
                    w = findall(vec(W) .== 1)
                    if isempty(w) || w[1] == 0
                        continue
                    end
                    behavior[i, "fieldIndex"] = w[1]
                end
                D[!, "fieldIndex"] = behavior[D[!,"nearest"], "fieldIndex"]
            end

            return D
        end

        """
        probability

        grabs the probability of the models given the data at each time
        Σ probabilityₜ(dataₜ)
        """
        function probability(data::DataFrame, fieldmodels::Dict)
            out =  Dict(name=>probability(data, vec(model), name) 
                        for (name,model) in fieldmodels)
            return out
        end
        function probability(data::DataFrame, fieldmodels::AbstractArray,
                units::AbstractVector)
            fieldIndex = data[!,:fieldIndex]
            probs = []
            for (model, unit) in zip(model, units)
                d = d[!, string(unit)]
                push!(probs, probability(d, model, fieldIndex))
            end
            return probs
        end
        function probability(data::DataFrame, fieldmodels::AbstractArray,
                name::NamedTuple)
            probability(data, fieldmodels, name.unit)
        end
        function probability(data::DataFrame, fieldmodel::AbstractVector,
                name::Int64)
            if string(name) in names(data)
                out = probability(data[!, string(name)], fieldmodel, data[!,:fieldIndex])
            else
                #println("Name=$name not in $(names(data))")
                out = nothing
            end
            return out
        end
        function probability(data::DataFrame, fieldmodels::AbstractArray,
                name::Int64)
            probability(data, fieldmodels[!, name], name)
        end
        function probability(counts::AbstractVector,
                fieldmodel::AbstractVector, fieldIndex::AbstractVector)
            goodInds = fieldIndex .> 0
            counts, fieldIndex = counts[goodInds], fieldIndex[goodInds]
            λ, k = fieldmodel[fieldIndex], counts
            return poisson.(k, λ)
        end
        function poisson(k::Real, λ::Real)
            out = (exp(-λ) * λ^k)/factorial(k)
            return out
        end

        """
         poisson_model : Creates a poisson_model from a field

        P⃗{x; λ}, where the model is valid ⟺   x is a count
        drawn from a single timebin Δt

        returns a function of field dimensions (to assign
        likelihood of each point)

        if grid == nothing
            then we assume users are passing in indices into our field distribution
        else
            users are passing a point, where we have to lookup the closted match
            in the grid
        end

        # modes
        K (give number of spikes), lambda is a Scalar
        lookup, lookup a behavior value of a spike, lambda is an Array
        lookup * K, K is a scalar and lambda is an Array

        """
        function poisson(k::Vector, λ::Vector, data::Tuple=nothing,
                grid::Tuple=nothing; lookup_K=true)
            if grid isa Vector && grid[1] isa Vector
                grid = Tuple(grid)
            end
            𝐊 = grid_lookup_func(grid)
            out = (exp.(-λ) .* λ.^(k * 𝐊(data)))./factorial.(k * 𝐊(data))
            return out
        end

        """
        grid_lookup_func

        creates a function that looks up field activity within a grid to
        and returns a count on that grid indicating number of events found
        per state.

        in essence, 

        it should take a point 𝐱 ∈ ℝᵈ and it will lookup it's location in
        the collection of edges 𝐆 ∈ ℝᵐ×ᵈ

        approach, we either:
        - fit() a histogram model from statsbase
        - nanstats
        """
        function grid_lookup_func(grid)
            lookup(data) = fit(Histogram, data, grid).weights
        end

        
    end
    export model

    module info
        using Entropies: Probabilities
        skipnan(x) = Iterators.filter(!isnan, x)

        """
        `spatial_information`

        computes the `spatial_information` 

        see yartsev dotson 2021 supp
        """
        function spatial_information(fields::Union{Dict,AbstractArray}, 
                behProb::Union{AbstractArray, Probabilities})
            if behfield isa AbstractArray
                behProb = Probabilities(collect(skipnan(vec(behfield))))
            end
            I = copy(fields)
            if fields isa Dict
                for i in keys(fields)
                    I[i] = spatial_information(fields)
                end
            else
                fields = collect(skipnan(vec(fields)))
                R = fields ./ mean(fields)
                I = sum( behProb .* R .* log2.(R) )
            end
            return I
        end

        function field_shift()
        end

        function field_shifts()
        end
    end
    export info


end

