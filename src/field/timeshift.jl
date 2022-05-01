module timeshift

    using DataFrames
    import ..field
    export get_field_shift
    export field
    using ThreadSafeDicts
    using DataStructures
    using DataFramesMeta
    using ProgressMeter
    using Distributions
    include("../table.jl")
    include("../shuffle.jl")
    import .table: to_dataframe
    using Infiltrator
    using Distributed
    using Dagger
    using Serialization
    using DrWatson

    export get_field_shift_shuffles, get_field_shift
    export to_dataframe, info_to_dataframe
    export fetch_best_fields
    export fetch_best_fields
    export shuffle_correct

    # -------------------- SHIFTING TYPES ---------------------------
    shift_func(data::DataFrame, shift::Real) = 
             transform(data, :time => (t->t.+shift) =>:time, copycols=false)
    const σ = shift_func

    # -------------------- SHIFTED Receptive Fields --------------------------
    function get_field_shift(beh::DataFrame, data::DataFrame, shift::Real; kws...)
        print("Using single")
        fieldobj = field.get_fields(σ(beh, shift), data; kws...)
    end

    function get_field_shift(beh::DataFrame, data::DataFrame,
            shifts::Union{StepRangeLen,Vector{T}} where T <: Real; 
            multithread::Bool=true,
            postfunc::Union{Function,Nothing}=nothing,
            as::Type=OrderedDict,
            multi::Union{Bool, Symbol}=true,
            safe_dict::AbstractDict=ThreadSafeDict(),
            kws...)

        kws = (;dokde=false, kws...)

        if multi isa Bool
            multi = multi ? :thread : :single
        end
        @info "Starting multi=$multi"
        msg = "$multi timeshift-shuffles"

        p = Progress(length(shifts), desc="field shift calculations")
        if multi == :thread
            Threads.@threads for shift in shifts
                if shift ∈ keys(safe_dict)
                    continue
                end
                result = field.get_fields(σ(beh, shift), data; kws...)
                if postfunc != nothing
                    result = postfunc(result)
                end
                push!(safe_dict, shift=>result)
                next!(p)
            end
        elseif multi == :single
            @showprogress for shift in shifts
                if shift ∈ keys(safe_dict)
                    continue
                end
                result = field.get_fields(σ(beh,shift), data; kws...)
                if postfunc != nothing
                    result = postfunc(result)
                end
                push!(safe_dict, shift=>result)
                next!(p)
            end
        elseif multi == :distributed
            @error "Not implemented"
        end
        safe_dict = Dict(safe_dict...)
        out = as(key=>pop!(safe_dict, key) for key in sort([keys(safe_dict)...]))
        return out
    end

    function get_field_shift_shuffles(beh::DataFrame, data::DataFrame,
            shifts::Union{StepRangeLen,Vector{T}} where T <: Real; 
            multi::Union{Symbol,Bool}=true,
            nShuffle::Int=100, 
            shuffle_pos::Union{Tuple,NamedTuple}=(;),
            shuffle_kws::NamedTuple=(;),
            shuffle_func=shuffle.by,
            postfunc::Union{Function,Nothing}=nothing,
            safe_dict::AbstractDict=ThreadSafeDict(),
            exfiltrateAfter::Real=Inf,
            kws...)::AbstractDict

        kws = (;dokde=false, kws...)

        sets = collect( Iterators.enumerate(Iterators.product(1:nShuffle,
                                                              shifts)))

        # Generate the distribution functoin for all shuffles (expensive to do for each)
        if isempty(shuffle_pos) || !(shuffle_pos[1] isa Distribution)
            shuffle_kws = (; shuffledist_df=beh, shuffle_kws...)
            distribution = shuffle._create_distribution(data, 
                                                        shuffle_pos...;
                                                        shuffle_kws...)
            shuffle_pos = (;shuffle_pos..., distribution)
            @info "distribution=$distribution"
        end
        
        if multi isa Bool
            multi = multi ? :thread : :single
        end

        @info "Starting multi=$multi"
        msg = "$multi timeshift-shuffles"
        skips = 0
        P = Progress(length(sets), desc=msg)
        if multi == :thread
            Threads.@threads for (i, (shuffle,shift)) in sets
                if (;shift,shuffle) ∈ keys(safe_dict)
                    skips+=1
                    continue
                end
                data = shuffle_func(data, shuffle_pos...; shuffle_kws...)
                result = field.get_fields(σ(beh,shift), data; kws...)
                if postfunc != nothing
                    result = postfunc(result)
                end
                push!(safe_dict, (;shift,shuffle)=>result)
                next!(P)
                if mod(i-skips, exfiltrateAfter) == 0
                    @info "chechpoint->exfiltrated"
                    @exfiltrate
                end
            end
        elseif multi == :distributed
            @assert nprocs() > 1 msg
            @showprogress 0.1 msg for (i,(shuffle,shift)) in sets
                data = Dagger.spawn(shuffle_func, data, shuffle_pos...; shuffle_kws...)
                @debug "Dagger 1"
                result = Dagger.spawn(field.get_fields, σ(beh,shift), data; kws...)
                @debug "Dagger 2"
                if postfunc != nothing
                    result = Dagger.@spawn postfunc(result)
                    @debug "Dagger 3"
                end
                push!(safe_dict, (;shift,shuffle)=>result)
            end
        elseif multi == :single
            for (i,(shuffle,shift)) in sets
                if (;shift,shuffle) ∈ keys(safe_dict)
                    skips+=1
                    next!(P)
                    continue
                else
                    #@info (; i,shift,shuffle)
                end
                data = shuffle_func(data, shuffle_pos...; shuffle_kws...)
                result = field.get_fields(σ(beh,shift), data; kws...)
                if postfunc != nothing
                    result = postfunc(result)
                end
                push!(safe_dict, (;shift,shuffle)=>result)
                if mod(i-skips, exfiltrateAfter) == 0
                    @info "chechpoint->exfiltrated"
                    @exfiltrate
                end
                next!(P)
            end
        else
            throw(ArgumentError("Unrecognized argument multi=$multi"))
        end
        safe_dict = Dict(safe_dict...)
        out = OrderedDict(key=>pop!(safe_dict, key) 
                          for key in sort([keys(safe_dict)...]))
        return out
    end


    function to_dataframe(shifts::AbstractDict{<:Real, <:Any}; kws...) 
        table.to_dataframe(shifts, key_name="shift", kws...)
    end
    function to_dataframe(
            shifts::AbstractDict{<:Union{NamedTuple,AbstractDict}, <:Any};
            kws...)
        table.to_dataframe(shifts, kws...)
    end

    function info_to_dataframe(shifts::AbstractDict{<:Real,<:Any};
            kws...)::DataFrame
        table.to_dataframe(shifts, key_name="shift", name="info", kws...)
    end
    function info_to_dataframe(
            shifts::AbstractDict{<:Union{NamedTuple,AbstractDict},<:Any};
            kws...)::DataFrame
        table.to_dataframe(shifts, name="info", kws...)
    end

    function fetch_best_fields(fieldInfo::DataFrame, pos...; kws...)
        beh, data = pos
        X = Dict()
        for neuron in fieldInfo.units
            D = @subset(data, :unit.==neuron)
            x = get_fields(σ(beh,fieldInfo.bestTau), D; kws...)
            push!(X,x)
        end
        return X
    end

    function shuffle_correct(main, shuffle)
    end

    function saveshifts(main, shuffle=nothing; metric=nothing, shifts=nothing,
                       fieldkws, tag="", kws...)

        parent_folder = datadir("exp_pro", "timeshift")
        if !(isdir(parent_folder))
            mkdir(parent_folder)
        end

        tag = isempty(tag) ? tag : "_$tag"
        if shifts != nothing
            start, stop = round(shifts[begin],digits=3),
                          round(shifts[end],  digits=3)
            N = length(shifts)
            shifts = "_Nstartstop=($N,$start:$stop)"
        else
            @error "No shifts name provided"
        end
        if metric == nothing
            @warn "No metric name provided, assuming metric='field'"
            metric = "field"
        end
        jf(x) = join(x,'-')
        props   = "props=$(jf(fieldkws.props))"
        splitby = "_splitby=$(jf(fieldkws.splitby))"

        name = joinpath(parent_folder, "$props$splitby$shifts$tag.serial")

        D = Dict(:main     => main,
                 :shuffle  => shuffle,
                 :shifts   => shifts,
                 :fieldkws => fieldkws)
        
        serialize(name, D)
        
    end


end
