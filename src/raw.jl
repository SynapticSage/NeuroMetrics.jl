module raw

    using DrWatson
    using Debugger
    using Plots
    using DataFrames
    using NetCDF
    using CSV
    using Statistics
    using Dates
    using Printf
    include(srcdir("utils", "SearchSortedNearest.jl", "src",
                       "SearchSortedNearest.jl"))
    animal_dayfactor = Dict("RY16"=>33, "RY22"=>0)
    export animal_dayfactor

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
            load_order = sort(data_source, by=source->source=="behavior",
                              rev=true)
        elseif as == "tuple"
            if "behavior" in data_source
                load_order = ("behavior", setdiff(data_source, ["behavior"])...)
            end
        else
            throw(ArgumentError("as cannot be $as"))
        end
        
        # Determine a time normalizing function
        if "behavior" ∈ data_source
            normalizing_time(data) = minimum(data["behavior"].time);
        else
            normalizing_time(data) = 0;
        end
        normalize(data, time) = (time .- normalizing_time(data))./60;

        # Load each data source
        data = Dict{String,Any}()
        for source ∈ load_order
            data[source] = raw.load_functions[source](args...)
            if "time" ∈ names(data[source])
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


    function load_behavior(animal::String, day::Int)
        function typeFunc(type, name)
            if occursin("Vec", string(name))
                type = ComplexF32;
            elseif name == "time"
                type = Float32;
            else
                type = nothing;
            end
            return type
        end
        typemap = Dict(Int64=>Int16);
        behCSV = datadir("exp_raw", "visualize_raw_neural",
                         "$(animal)_$(day)_beh.csv")
        beh = CSV.read(behCSV, DataFrame,
                 strict=false,
                 missingstring=["NaNNaNi", "NaNNaNi,", ""],
                 types=typeFunc,
                 typemap=typemap)
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

    function load_spikes(animal::String, day::Int; beh=Nothing)
        rawSpikingCSV = DrWatson.datadir("exp_raw",
                                         "visualize_raw_neural",
                                         "$(animal)_$(day)_labeled_spiking.csv")
        raster = CSV.read(rawSpikingCSV, DataFrame,
                 strict=false,
                 missingstring=["NaN", "", "NaNNaNi"])


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

    function load_cells(animal::String, day::Int)
        csvFile = DrWatson.datadir("exp_raw",
                                         "visualize_raw_neural",
                                         "$(animal)_$(day)_cell.csv")
        cells = CSV.read(csvFile, DataFrame,
                 strict=false,
                 missingstring=["NaN", "", "NaNNaNi"])
        return cells
    end

    function load_ripples(animal, day)
        typemap = Dict(Int64=>Int16);
        csvFile = DrWatson.datadir("exp_raw",
                                   "visualize_raw_neural",
                                   "$(animal)_$(day)_ripple.csv")
        ripples = CSV.read(csvFile, DataFrame, strict=false,
                 typemap=typemap,
                 missingstring=["NaN", "NaNNaNi", "NaNNaNi,", ""])
    end

    function load_lfp(animal, day; vars=nothing)
        netcdf = DrWatson.datadir("exp_raw",
                                         "visualize_raw_neural",
                                         "$(animal)_$(day)_rhythm.nc")
        v = NetCDF.open(netcdf)
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

    function load_task(animal, day)
        typemap = Dict(Int64=>Int16);
        csvFile = DrWatson.datadir("exp_raw",
                                   "visualize_raw_neural",
                                   "$(animal)_$(day)_task.csv")
        task = CSV.read(csvFile, DataFrame, strict=false,
                 typemap=typemap,
                 missingstring=["NaN", "NaNNaNi", "NaNNaNi,", ""])
        transform!(task, :x => pxtocm => :x, :y => pxtocm => :y)
        return task
    end

    function decodedir(method="sortedspike", transition="empirical", binstate="notbinned", n_split=4, downsamp=1, speedup=1.0)
        base = "/Volumes/FastData/decode/"
        paramfolder = "$method.$transition.$binstate.n_split=$n_split.downsamp=$downsamp.speedup=$(@sprintf("%1.1f", speedup))"
        return joinpath(base, paramfolder)
    end
    function decodepath(animal="RY16", day=36, epoch=7; type="test", split=1, kws...)
        dir = decodedir(kws...)
        file = "$(animal)decode$(@sprintf("%02d",day))-$(@sprintf("%02d", epoch))split=$split.$type.nc"
        joinpath(dir, file)
    end

    function load_decode(filename)
        v = NetCDF.open(filename)
        if "Var1" in keys(v.vars)
            v.vars["time"] = v.vars["Var1"]
            pop!(v.vars, "Var1")
        end
        lfp = Dict(var => Array(v.vars[var]) 
                   for var in keys(v.vars))
        return lfp
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
        return tetrodes[!, Not(:cell)]
    end

    function load_pathtable(animal, day)
        f=CSV.read(datadir("paths.csv"), DataFrame)
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
    
    module behavior
        function add_next_target(beh)

        end
        function add_previous_target(beh)
        end
        additions = Dict("previous_target" => add_previous_target,
                         "next_target"=> add_next_target)
    end
    export behavior

    module task
        function well_locations(task; day=nothing, epoch=nothing)
        end
    end
    export task

    module lfp
        function annotate_cycles(lfp; phase_col="phase")
            phase = lfp[!, phase_col]
            Δₚ = [0; diff(phase)]
            change_points = Δₚ .< 0
            cycle_labels = accumulate(+, change_points)
        end
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
        (downsamp)
    end

    function filterTables(data::DataFrame...;
            filters::Union{Nothing,Dict}=Dict(),
            lookupcols=nothing, lookupon="time")::Vector{DataFrame}

        if filters == nothing
            filters = Dict()
        end

        data = [data...];
        # Get our columns into target
        if lookupcols != nothing
            for col ∈ lookupcols
                source = col[1].source
                target = col[1].target
                columns_to_transfer   = col[2]
                if columns_to_transfer == All()
                    continue
                end

                match_on_source = data[source][:, lookupon]
                match_on_target = data[target][:, lookupon]
                match_on_target = convert.(Float64, match_on_target)
                match_on_source = convert.(Float64, match_on_source)
                match_on_source = (match_on_source,)
                indices_of_source_samples_in_target = SearchSortedNearest.searchsortednearest.(match_on_source, match_on_target)
                for (i, item) ∈ Iterators.enumerate(columns_to_transfer)
                    data[target][!, item] = 
                        data[source][indices_of_source_samples_in_target, item]
                end
            end
        end

        # Filtration
        println("→ → → → → → → → → → → → ")
        println("Filtration")
        println("→ → → → → → → → → → → → ")
        for filt ∈ filters
            for i ∈ 1:length(data)
                #print(filt)
                filter_cols, filter_function = filt
                @assert !(filter_cols isa Bool)
                @assert !(filter_function isa Bool)
                @assert !(filter_function isa Vector{Bool})
                @assert !(filter_cols isa Vector{Bool})
                inds = filter_function(data[i][!, filter_cols]);
                percent = mean(inds)*100
                println("data_$i filtration: $percent percent pass filter")
                x = data[i][findall(inds), :];
                data[i] = x;
            end
        end
        println("← ← ← ← ← ← ← ← ← ← ← ← ")

        return data

    end

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
                dayfactor = animal_dayfactor[animal]
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
            stream = VideoIO.open(videopath)
            vid = VideoIO.openvideo(stream)
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
            df = CSV.read(get_path(pos...;kws...), DataFrame, header=3, skipto=4)
            transform!(df, :coords=>(x->x*(1/30))=>:time, 
                          [:x,:x_1]=>((a,b)->(a.+b)./2)=>:X,
                          [:y,:y_1]=>((a,b)->(a.+b)./2)=>:Y,
                   )

            transform!(df, [:X, :Y]=>((X,Y)->sqrt.(X.^2 .+ Y.^2))=>:vec)
            transform!(df, :vec => (x->[0;diff(x)]) => :delta)
        end
    end
    export dlc

end
