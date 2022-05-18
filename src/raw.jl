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

    include("./raw/behavior.jl")
    include("./raw/lfp.jl")
    include("./raw/task.jl")
    include("./raw/decode.jl")
    include("./raw/celltet.jl")
    include("./raw/spikes.jl")
    include("./raw/ripples.jl")

    load_functions = Dict(
        "behavior" => load_behavior,
        "spikes"   => load_spikes,
        "lfp"      => load_lfp,
        "cells"    => load_cells,
        "tetrode"  => load_tetrode,
        "decode"   => load_decode,
        "task"     => load_task,
        "ripples"  => load_ripples)

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
    
    # ----------
    # OPERATIONS
    # ----------

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
