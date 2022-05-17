module raw
    
    #Imports
    using DrWatson
    using DataFrames
    using NetCDF
    using CSV
    using Statistics
    using Dates
    using Printf
    using ProgressMeter
    using Glob
    using Infiltrator
    include("utils.jl")
    import .utils
    include("table.jl")
    import .table
    using .table: get_periods
    const findnearest = utils.searchsortednearest
    using Infiltrator
    __revise_mode__ = :evalassign

    # Module-wide settings
    animal_dayfactor = Dict("RY16"=>33, "RY22"=>0)
    csvkws=(; silencewarnings=true, buffer_in_memory=true, ntasks=1)
    pxtocm(x) = 0.1487 * x
    cmtopx(x) = x / 0.1487 

    """
        load(animal, day)

    Loads up the datasets we plan to use to plot out the raw raster data

    =====
    Input
    =====
    animal : String
        Name of the animal
    day : Int
        Integer of the day


    ======
    Output
    ======
    (raster, behavior)

    """
    function load(args...; as="tuple",
            data_source=["spikes", "behavior", "ripples", "cells"])

        # Establish a load order
        if as == "dict"
            load_order = sort(data_source, by=source->source=="behavior", rev=true)
        elseif as == "tuple"
            if "behavior" in data_source
                load_order = ("behavior", setdiff(data_source, ["behavior"])...)
            end
        else
            throw(ArgumentError("as cannot be $as"))
        end
        
        # Determine a time normalizing function
        if "behavior" ∈ data_source
            normalizing_time(data) = minimum(data["behavior"].time)
        else
            normalizing_time(data) = 0
        end

        normalize(data, time) = (time .- normalizing_time(data))./60;

        # Load each data source
        data = Dict{String,Any}()
        for source ∈ load_order
            data[source] = raw.load_functions[source](args...)
            if "time" ∈ names(data[source])
                @info "Centering $source and **converting to minutes**"
                data[source].time = normalize(data, data[source].time);
            end
        end

        # How do we return, dict or tuple?
        if as == "tuple"
            return (data[name] for name in data_source)
        else
            return data
        end
        
    end


    function load_spikes(animal::String, day::Int; beh=Nothing)
        rawSpikingCSV = DrWatson.datadir("exp_raw",
                                         "visualize_raw_neural",
                                         "$(animal)_$(day)_labeled_spiking.csv"
                                        )
        @info rawSpikingCSV
        raster = CSV.read(rawSpikingCSV, DataFrame;
                 strict=false, missingstring=["NaN", "", "NaNNaNi"],
                 csvkws...)

        # And let's add some brain area specific numbering
        groups = groupby(raster, "area");
        for g = 1:length(groups)
            unit = unique(groups[g].unit)
            areaunit = 1:length(unit);
            groups[g].areaunit = map(x -> Dict(unit .=> areaunit)[x],
                                     groups[g].unit)
        end
        raster = combine(groups, x->x)
    end

    function load_ripples(animal, day)
        typemap = Dict(Int64=>Int16);
        rippleFile = DrWatson.datadir("exp_raw",
                                   "visualize_raw_neural",
                                   "$(animal)_$(day)_ripple.csv")
        @info rippleFile
        ripples = CSV.read(rippleFile, DataFrame; strict=false,
                 typemap=typemap,
                 missingstring=["NaN", "NaNNaNi", "NaNNaNi,", ""],
                 csvkws...)
    end

    include("./raw/behavior.jl")
    include("./raw/lfp.jl")
    include("./raw/task.jl")
    include("./raw/decode.jl")
    include("./raw/cells.jl")


    function load_tetrode(animal,day)
        cells = load_cells(animal,day)
        groups = groupby(cells,"tetrode")
        tetrodes = DataFrame()
        for group = groups
            n_cells = size(group,1)
            row = DataFrame(group[1,:])
            row[!, :n_cells] .= n_cells;
            append!(tetrodes, row);
        end
        out = if "cell" in names(tetrodes)
                tetrodes[!, Not(:cell)]
            else
                tetrodes
            end
        return out
    end

    function load_pathtable(animal, day)
        f=CSV.read(datadir("paths.csv"), DataFrame; csvkws...)
        init = Dates.Time("00:00:00")
        function s(x)
            x = Second.(x .- init)
            x = [e.value for e in x]
        end
        transform!(f, 
                   :start=>(x->s(x))=>:start, 
                   :end=>(x->s(x))=>:end, 
                   [:end,:start]=>((e,s)->Dates.Second.(e.-s))=>:duration)
        return f
    end

    load_functions = Dict(
        "behavior" => load_behavior,
        "spikes"   => load_spikes,
        "lfp"      => load_lfp,
        "cells"    => load_cells,
        "tetrode"  => load_tetrode,
        "decode"   => load_decode,
        "task"     => load_task,
        "ripples"  => load_ripples)
    
    # ----------
    # OPERATIONS
    # ----------

    """
        downsample(animal, day)

    Downsamples the data

    Input
    =====
    raster : DataFrame
    beh : DataFrame


    Output
    ======
    (raster, behavior)

    """
    function downsample(x; dfactor=10)
        δ = dfactor;
        downsamp = x[begin:δ:end, :];
        downsamp
    end

    """
    function register(data::DataFrame...; transfer,
            on::String="time")::Vector{DataFrame} 
    """
    function register(data::DataFrame...; transfer,
            on::String="time")::Vector{DataFrame} 
        if data isa Tuple
            data = [data...];
        end
        # Get our columns into target
        @debug "→ → → → → → → → → → → → "
        @debug "Registration"
        @debug "→ → → → → → → → → → → → "
        for col ∈ transfer
            source = col[1].source
            target = col[1].target
            columns_to_transfer   = col[2]
            if columns_to_transfer == All()
                continue
            end
            @debug "columns=$columns_to_transfer from source->target on $on"

            #@infiltrate
            match_on_source = data[source][:, on]
            match_on_target = data[target][:, on]
            match_on_target = convert.(Float64, match_on_target)
            match_on_source = convert.(Float64, match_on_source)
            match_on_source = (match_on_source,)
            indices_of_source_samples_in_target = findnearest.(match_on_source, match_on_target)
            for (i, item) ∈ Iterators.enumerate(columns_to_transfer)
                data[target][!, item] =
                data[source][indices_of_source_samples_in_target, item]
            end
        end
        @debug "← ← ← ← ← ← ← ← ← ← ← ← "
        return data
    end
    """
    register

    function register(source::DataFrame, target::DataFrame; 
            transfer, on::String="time")::Vector{DataFrame}

    register columns in a `source` to a `target` dataframe `on` a certain
    column
    """
    function register(source::DataFrame, target::DataFrame; 
            transfer, on::String="time")::Union{Tuple, DataFrame}
        if transfer isa String
            transfer = [transfer]
        end
        if transfer isa Vector{String}
            addressing = (;source=1, target=2)
            transfer = ((addressing, transfer),) # create set of addressed transfer instructions
        end
        @debug "Got here"
        source, target, _ = register(source, target, DataFrame(); transfer=transfer, on=on)
        return source, target
    end
    function registerEventsToContinuous(events::DataFrame, target::DataFrame;
            transfer::Union{Vector{String},String}, on::String="time",
            targetEltype::Dict{String,<:Type}=Dict{String,Any}(),
            targetCast::Dict{String,<:Function}=Dict{String,Function}(),
            ifNonMissingAppend::Bool=false,
            eventStart::String="start",
            eventStop::String="stop")::DataFrame 

        # Get our columns into target
        @debug "→ → → → → → → → → → → → "
        @debug "Registration"
        @debug "→ → → → → → → → → → → → "
        
        for col ∈ transfer
            if col ∉ names(target)
                @debug "$col not in target"
                T = col ∈ keys(targetEltype) ? targetEltype[col] : eltype(events[!,col])
                @debug "col=$col => type=$T"
                @debug "creating new cool with type = $T"
                target[!,col] = Array{Union{Missing,T}}(missing, size(target,1))
                # Cast to user specified type?
                if col in keys(targetCast)
                    target[!,col] = targetCast[col](target[!,col])
                end
            end
        end

        match = target[!, on]
        p = Progress(size(events,1), desc="Registering")
        Threads.@threads for event in eachrow(events)
            start = event[eventStart]
            stop  = event[eventStop]
            indices_of_source_samples_in_target = (match .>= start) .&&
                                                  (match .< stop)
            #@debug "sum=$(sum(indices_of_source_samples_in_target))"
            for col ∈ transfer
                if col == All()
                    continue
                end
                specialStringBehavior = ifNonMissingAppend && 
                         eltype(target[!, col]) == Union{Missing,String} # has to be String, not InlineString
                #@debug "$col is eltype=$(eltype(target[!,col])) and specialStringBehavior=$specialStringBehavior"
                if specialStringBehavior
                    #@debug "special behavior"
                    missingVals    = indices_of_source_samples_in_target .&& ismissing.(target[!,col])
                    nonMissingVals = indices_of_source_samples_in_target .&& (!).(ismissing.(target[!,col]))
                    existing       = target[nonMissingVals, col]
                    target[nonMissingVals, col] .= existing .* event[col]
                    target[missingVals, col]    .= event[col]
                else
                    target[indices_of_source_samples_in_target, col] .= event[col]
                end
            end
            next!(p)
        end
        @debug "← ← ← ← ← ← ← ← ← ← ← ← "
        return target
    end
    """
    filter

    instructions to query/filter values in a set of dataframes
    """
    function filter(data::DataFrame...; filters::AbstractDict=Dict())::Vector{DataFrame}
        data = [data...]
        @debug "→ → → → → → → → → → → → "
        @debug "Filtration"
        @debug "→ → → → → → → → → → → → "
        for (cols, filt_for_cols) ∈ filters
            for i ∈ 1:length(data)
                @assert !(cols isa Bool)
                @assert !(filt_for_cols isa Bool)
                @assert !(filt_for_cols isa Vector{Bool})
                @assert !(cols isa Vector{Bool})
                if filt_for_cols isa Function
                    #println("Filter is a function")
                    inds = filt_for_cols(data[i][!, cols]);
                elseif typeof(filt_for_cols) <: Vector
                    #println("Filter is a set of functions")
                    inds = accumulate(.&, [ff(data[i][!, cols]) for ff
                                           in filt_for_cols])[end]
                else
                    @debug "cols = $cols"
                    @debug "Typeof(cols) = $(typeof(cols))"
                    @debug "filt_for_cols = $filt_for_cols"
                    throw(TypeError(:filt_for_cols, "",Vector,typeof(filt_for_cols)))
                end
                percent = mean(inds)*100
                @debug "data_$i filtration: $percent percent pass filter on $cols"
                x = data[i][findall(inds), :];
                data[i] = x;
            end
        end
        for i in 1:length(data)
            @debug "final_size(data($i)) = $(size(data[i]))"
        end
        @debug "← ← ← ← ← ← ← ← ← ← ← ← "
        return data
    end

    """
    `filterAndRegister`

    combination of  `raw.filter()` and `raw.register()` steps
    """
    function filterAndRegister(data::DataFrame...;
            filters::Union{Nothing,AbstractDict}=nothing, transfer=nothing,
            on="time")::Vector{DataFrame}
        if transfer != nothing
            data = register(data...; transfer=transfer, on=on)
            #utils.piso(data)
        end
        if filters != nothing
            data = filter(data...; filters=filters)
        end
        return data    
    end

    """
    `filterTables`

    alias for `filterAndRegister`
    """
    filterTables(data::DataFrame...; filters::Union{Nothing,AbstractDict}=nothing,
                 lookupcols=nothing, lookupon="time") =
    filterAndRegister(data...;transfer=lookupcols, on=lookupon, filters=filters)

    module video
        using Glob, Printf
        using VideoIO
        using ..raw 
        function frameattime(vid, time; cropx=[], cropy=[])
            if time == 0
                vid = seekstart(vid)
            else
                currtime = gettime(vid)
                vid = seek(vid, time - currtime)
            end
            seek(vid, time)
            img = read(vid)'
            if (length(cropx) & length(cropy)) > 0
                cropx = Int.(round.(pxtocm.(cropx)))
                cropy = Int.(round.(pxtocm.(cropy)))
                img=img[cropx[1] : cropx[2], 
                    cropy[1] : cropy[2]]
            end
            #img = img[:, end:-1:begin]
            img
        end
        function get_path(animal, day, epoch; dayfactor=0, 
                guessdayfactor=true,
                source="deeplabcut")
            if guessdayfactor
                dayfactor = raw.animal_dayfactor[animal]
            end
            day += dayfactor
            if source == "deeplabcut"
                folder_path = "/Volumes/Colliculus/deeplabcut/" * 
                             "goalmaze_tape-Ryan-2020-05-28/videos/"
            end
            possible_files = 
            glob("$(animal)_$(@sprintf("%02d",day))_$(@sprintf("%02d",epoch))_*.mp4",
                     folder_path)
            videopath = possible_files[1]
            return videopath
        end
        function load(pos...; kws...)
            videopath = get_path(pos...; kws...)
            stream    = VideoIO.open(videopath)
            vid       = VideoIO.openvideo(stream)
        end
    end
    export video

    module dlc
        using Glob, Printf
        using CSV, DataFrames
        using ..raw 
        function get_path(animal, day, epoch; dayfactor=0, 
                guessdayfactor=true, filtered=false, source="deeplabcut")
            if guessdayfactor
                dayfactor = raw.animal_dayfactor[animal]
            end
            day += dayfactor
            if source == "deeplabcut"
                folder_path = "/Volumes/Colliculus/deeplabcut/" * 
                             "goalmaze_tape-Ryan-2020-05-28/videos/"
            end
            possible_files = 
                glob("$(animal)_$(@sprintf("%02d",day))_$(@sprintf("%02d",epoch))_*.csv",
                     folder_path)
            #print(possible_files)
            if filtered
                possible_files = [file for file in possible_files if
                                  occursin(file,"filtered")]
            else
                possible_files = [file for file in possible_files if
                                  !(occursin(file,"filtered"))]
            end
            videopath = possible_files[1]
            return videopath
        end
        function load(pos...; kws...)
            df = CSV.read(get_path(pos...;kws...), DataFrame; header=3, skipto=4,
                         csvkws...)
            transform!(df, :coords=>(x->x*(1/30))=>:time, 
                          [:x,:x_1]=>((a,b)->(a.+b)./2)=>:X,
                          [:y,:y_1]=>((a,b)->(a.+b)./2)=>:Y,
                   )

            transform!(df, [:X, :Y]=>((X,Y)->sqrt.(X.^2 .+ Y.^2))=>:vec)
            transform!(df, :vec => (x->[0;diff(x)]) => :delta)
        end
    end
    export dlc

    function normalize_time(data::Union{DataFrame, Dict}...;
            timefields::Dict=Dict(), factor=1, warnifnotime::Bool=false)
        data = [data...]
        if data[1] isa DataFrame
            tₘ = minimum(data[1].time)
        else
            tₘ = minimum(data[1]["time"])
        end
        for source ∈ 1:length(data)
            @debug "source = $source"
            if !(source in keys(timefields))
                tfs = ["time"]
                @debug "tfs = $tfs"
            else
                tfs = timefields[source]
                @debug "tfs = $tfs"
            end
            for tf in tfs
                if data[source] isa DataFrame && (tf ∈ names(data[source]))
                    data[source][!,tf] = (data[source][!,tf] .- tₘ)*factor
                elseif !(data[source] isa DataFrame) && (tf ∈ keys(data[source]))
                    data[source][tf] = (data[source][tf] .- tₘ)*factor
                else
                    if warnifnotime
                        @warn "No time_field=$tf for $source"
                    else
                        @error "No time_field=$tf for $source"
                    end
                end
            end
        end
        return data
    end

    function keep_overlapping_times(data::Union{DataFrame, Dict, Vector}...; 
            tf="time", returninds::Vector=[])

        data=[data...]
        
        # Get time extrema
        sz = []
        for d in data
            if d isa DataFrame && (tf ∈ names(d))
                push!(sz, extrema(d[!,tf]))
            elseif !(d isa DataFrame) && (tf ∈ keys(d))
                push!(sz, extrema(d[tf]))
            elseif d isa Vector
                push!(sz, extrema(d)) 
            end
        end

        # Compute overall extrema
        sz = cat([[s...] for s in sz]...; dims=2)'
        minmax = [maximum(sz[:,1]), minimum(sz[:,2])]


        # Contrain each object to live in the overall extrema
        ind_constrain(x) = (x .>= minmax[1]) .&& (x .< minmax[2])
        for (i, d) in zip(1:length(data), data)
            if d isa DataFrame && (tf ∈ names(d))
                inds = ind_constrain(d[!,tf])
            elseif !(d isa DataFrame) && (tf ∈ keys(d))
                inds = ind_constrain(d[tf]) 
            elseif d isa Vector
                inds = ind_constrain(d) 
            end

            
            if i in returninds
                data[i] = inds
            else
                if d isa DataFrame && (tf ∈ names(d))
                    data[i] = d[inds, :]
                elseif !(d isa DataFrame) && (tf ∈ keys(d))
                    time_length = length(d[tf])
                    dims = Dict(key=>findfirst(size(value).==time_length)
                                for (key,value) in d)
                    for (k,v) in d
                        I = Vector{Any}([Colon() for i in 1:ndims(v)])
                        if !(isempty(dims[k]))
                            I[dims[k]] = inds
                        end
                        data[i][k] = getindex(v, I...)
                    end
                elseif d isa Vector
                    data[i] = d[inds]
                end
            end
        end

        return data
    end

end
