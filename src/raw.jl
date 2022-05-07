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
    const findnearest = utils.searchsortednearest

    # Module-wide settings
    animal_dayfactor = Dict("RY16"=>33, "RY22"=>0)
    csvkws=(; silencewarnings=true, buffer_in_memory=true, ntasks=1)

    """
        load(animal, day)

    Loads up the datasets we plan to use to plot out the raw raster data

    Input
    =====
    animal : String
        Name of the animal
    day : Int
        Integer of the day


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


    function behaviorpath(animal::String, day::Int, tag::String="")
        tag  = length(tag) == 0 ? tag : "_$tag"
        path = datadir("exp_raw", "visualize_raw_neural",
                         "$(animal)_$(day)_beh$tag")
        @debug "path=$path"
        if occursin("*",tag)
            path = glob(basename(path), dirname(path))
        end
        @debug "path=$path"
        return path
    end
    function load_behavior(animal::String, day::Int, tag::String="")
        function typeFunc(type, name)
            if occursin("Vec", string(name))
                type = ComplexF32;
            elseif name == "time" type = Float32;
            else
                type = nothing;
            end
            return type
        end
        typemap = Dict(Int64=>Int16);
        @debug "animal=$animal, day=$day, tag=$tag"
        behCSV = behaviorpath(animal, day, tag) .* ".csv"
        @info "behCSV=>$behCSV"
        readFile(file) = CSV.read(file, DataFrame;
                                  strict=false, 
                                  missingstring=["NaNNaNi", "NaNNaNi,", ""], 
                                  types=typeFunc, typemap=typemap, csvkws...)
        if behCSV isa Vector
            beh = [readFile(file) for file in behCSV]
            beh = hcat(beh...)
        else
            beh = readFile(behCSV)
        end
        if beh.time isa Vector{String}
            beh.time = parse.(Float32, beh.time);
        end
        for col in names(beh)
            if occursin("current", col) || occursin("egoVec", col)
                replace!(beh[!,col], missing=>NaN)
            end
        end
        @assert ("x" ∈ names(beh)) "Fuck"
        return beh
    end
    function save_behavior(animal::String, day::Int)
    end

    function load_spikes(animal::String, day::Int; beh=Nothing)
        rawSpikingCSV = DrWatson.datadir("exp_raw",
                                         "visualize_raw_neural",
                                         "$(animal)_$(day)_labeled_spiking.csv"
                                        )
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

    function cellpath(animal::String, day::Int, tag::String=""; kws...)
        if tag != "" && tag != "*"
            tag = "_$tag"
        end
        csvFile = DrWatson.datadir("exp_raw", "visualize_raw_neural",
                                   "$(animal)_$(day)_cell$tag.csv")
    end
    function load_cells(pos...; kws...)
        path = cellpath(pos...; kws...)
        if occursin("*", path)
            base, dir = basename(path), dirname(path)
            @debug "base=$base, dir=$dir"
            paths = glob(base, dir)
        else
            paths = [path]
        end

        cells = DataFrame()
        @showprogress 0.1 "loading cell files" for path in paths
            cell = CSV.read(path, DataFrame;
                     strict=false, missingstring=["NaN", "", "NaNNaNi"],
                     csvkws...)
            cells = isempty(cells) ? cell : outerjoin(cells, cell, on=:unit, makeunique=true)
            table.clean_duplicate_cols(cells)
        end
        return cells
    end
    function save_cells(cells::DataFrame, pos...; merge_if_exist::Bool=true, kws...)
        #kws = (;kws...) # satellite default is true
        csvFile = cellpath(pos...; kws...)
        #if merge_if_exist && isfile(csvFile)
        #    prevcells = load_cells(pos...; kws...)
        #    cells = outerjoin(cells, prevcells, makeunique=true)
        #    # TODO function to accept left or right dups
        #end
        println("Saving cell data at $csvFile")
        cells |> CSV.write(csvFile)
    end

    function load_ripples(animal, day)
        typemap = Dict(Int64=>Int16);
        csvFile = DrWatson.datadir("exp_raw",
                                   "visualize_raw_neural",
                                   "$(animal)_$(day)_ripple.csv")
        ripples = CSV.read(csvFile, DataFrame; strict=false,
                 typemap=typemap,
                 missingstring=["NaN", "NaNNaNi", "NaNNaNi,", ""],
                 csvkws...)
    end

    function lfppath(animal::String, day::Int)
        netcdf = DrWatson.datadir("exp_raw", "visualize_raw_neural",
                                         "$(animal)_$(day)_rhythm.nc")
    end
    function load_lfp(pos...; vars=nothing)
        v = NetCDF.open(lfppath(pos...))
        if "Var1" in keys(v.vars)
            v.vars["time"] = v.vars["Var1"]
            pop!(v.vars, "Var1")
        end
        keyset = keys(v.vars)
        if vars != nothing
            throw(InvalidStateException("Not impllemented yet"))
        end
        lfp = Dict(var => Array(v.vars[var]) 
                   for var in keyset)
        lfp = DataFrame(Dict(var => vec(lfp[var]) 
                             for var in keyset))
        return lfp
    end
    function cyclepath(animal::String, day::Int, tetrode::Union{String,Int})
        csv = DrWatson.datadir("exp_raw", "visualize_raw_neural",
                                  "$(animal)_$(day)_tet=$(tetrode)_cycles.csv")
    end
    function save_cycles(cycles, pos...)
        cycles |> CSV.write(cyclepath(pos...))
    end
    function load_cycles(pos...)
        cycles = CSV.read(cyclepath(pos...), DataFrame; strict=false,
                 missingstring=["NaN", "NaNNaNi", "NaNNaNi,", ""], csvkws...)
    end

    function load_task(animal::String, day::Int)
        typemap = Dict(Int64=>Int16);
        csvFile = DrWatson.datadir("exp_raw",
                                   "visualize_raw_neural",
                                   "$(animal)_$(day)_task.csv")
        task = CSV.read(csvFile, DataFrame; strict=false, typemap=typemap,
                 missingstring=["NaN", "NaNNaNi", "NaNNaNi,", ""], csvkws...)
        transform!(task, :x => pxtocm => :x, :y => pxtocm => :y)
        return task
    end

    function decodedir(;method::String="sortedspike", transition="empirical",
            binstate="notbinned", n_split=4, downsamp=1, speedup=1.0)
        base = "/Volumes/FastData/decode/"
        paramfolder = "$method.$transition.$binstate.n_split=$n_split.downsamp=$downsamp.speedup=$(@sprintf("%1.1f", speedup))"
        return joinpath(base, paramfolder)
    end
    function decodepath(animal="RY16", day=36, epoch=7; type="test", split=1, kws...)
        dir = decodedir(;kws...)
        file = "$(animal)decode$(@sprintf("%02d",day))-$(@sprintf("%02d", epoch))split=$split.$type.nc"
        fullpath = joinpath(dir, file)
        if !(isfile(fullpath))
            @warn "file=$fullpath does not exist"
        end
        return fullpath
    end

    function load_decode(filename)
        v = NetCDF.open(filename)
        if "Var1" in keys(v.vars)
            v.vars["time"] = v.vars["Var1"]
            pop!(v.vars, "Var1")
        end
        decode = Dict(var => Array(v.vars[var]) 
                   for var in keys(v.vars))
        return decode
    end

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
    

    module task
        function well_locations(task; day=nothing, epoch=nothing)
        end
    end
    export task

    module lfp
        using DataFrames
        using Statistics
        using DirectionalStatistics
        using ImageFiltering
        include("table.jl")
        function phase_to_radians(phase)
            phase = Float32.(phase)
            phase = 2*π*(phase .- minimum(phase))./diff([extrema(phase)...]) .- π
            #phase = convert.(Float16, phase)
        end
        function annotate_cycles(lfp; phase_col="phase", method="peak-to-peak")
            phase = lfp[!, phase_col]
            lfp.phase = phase_to_radians(lfp[:,"phase"])
            println("Method=$method")
            if method == "resets"
                Δₚ = [0; diff(phase)]
                reset_points = UInt32.(Δₚ .< (-1.5*π))
                cycle_labels = accumulate(+, reset_points)
            elseif method == "peak-to-peak"
                step_size = median(diff(phase))
                Δₚ = [0; diff(phase)]
                #falling_zero_point = [(phase[1:end-1] .>=0) .& (phase[2:end] .<0) ; false]
                rising_zero_point = [(phase[2:end] .>=0) .& (phase[1:end-1] .<0) ; false]
                cycle_labels = accumulate(+, rising_zero_point)
                lfp[!,"phase"] = mod2pi.(lfp[!,"phase"])
            elseif method == "trough-to-trough"
                step_size = median(diff(phase))
                Δₚ = [0; diff(phase)]
                falling_zero_point = [(phase[1:end-1] .>=0) .& (phase[2:end] .<0) ; false]
                #rising_zero_point = [(phase[2:end] .>=0) .& (phase[1:end-1] .<0) ; false]
                cycle_labels = accumulate(+, falling_zero_point)
                lfp[!,"phase′"] = mod.(lfp[!,"phase"] .- pi, 2*pi)
            else
                throw(ArgumentError("Unrecognized method=$method"))
            end
            lfp[!,"cycle"] = cycle_labels
            return lfp
        end
        function mean_lfp(lfp; mean_fields=["phase","amp","raw"], func=Circular.median)
            lfp = groupby(lfp, :time)
            non_mean_fields = setdiff(names(lfp), mean_fields)
            lfp =combine(lfp, mean_fields.=>func.=>mean_fields, 
                         non_mean_fields.=>first.=>non_mean_fields)
            return lfp[!, Not(:tetrode)]
        end
        function gauss_lfp(lfp; fields=["phase"], gaussian=3)
            kernel = Kernel.gaussian((gaussian,))
            for field in fields
                lfp[!,field] = imfilter(lfp[!,field], kernel)
            end
            return lfp
        end
        """
        weighted_lfp

        creates an amplitude weighted average of fields across tetrodes
        """
        function weighted_lfp(lfp; mean_fields=["phase","amp","raw"],
                weighting="amp")
            lfp = groupby(lfp, :time)
            non_mean_fields = setdiff(names(lfp), mean_fields)
            new = DataFrame()
            @time for lf in lfp
                item = DataFrame(lf[1,:])
                for field in mean_fields
                    item[!, field] .= sum(lf[!, field] .* lf[!, weighting])/sum(lf.amp)
                end
                append!(new, item)
            end
            return new[!, Not(:tetrode)]
        end

        function unstack_tetrode(df; measure::Symbol=:phase)
            unstack(df[!, [:time, :tetrode, measure]], :tetrode, measure)
        end

        function get_cycle_table(lfp, pos...; kws...)
            @assert "cycle" in names(lfp)
            tab = table.get_periods(lfp, "cycle", :amp=>mean, pos...; kws...)
            return tab
        end
        getTet(L::DataFrame, T::Int) = filter(:tetrode=> t->t==T, L)

    end
    export lfp


    # ----------
    # OPERATIONS
    # ----------
    pxtocm(x) = 0.1487 * x
    cmtopx(x) = x / 0.1487 

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

            @infiltrate
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
