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
    export table
    using Infiltrator
    using Distributed
    using Dagger

    # -------------------- SHIFTING TYPES ---------------------------
    shift_func(data, shift) = transform(data, :time => (t->t.+shift) =>:time, copycols=false)
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
            kws...)

        if multi isa Bool
            multi = multi ? :thread : :single
        end
        @info "Starting multi=$multi"
        msg = "$multi timeshift-shuffles"

        safe_dict = ThreadSafeDict()
        p = Progress(length(shifts), desc="field shift calculations")
        if multi
            Threads.@threads for shift in shifts
                result = field.get_fields(σ(beh,shift), data; kws...)
                if postfunc != nothing
                    result = postfunc(result)
                end
                push!(safe_dict, shift=>result)
                next!(p)
            end
        else
            @showprogress for shift in shifts
                result = field.get_fields(σ(beh,shift), data; kws...)
                if postfunc != nothing
                    result = postfunc(result)
                end
                push!(safe_dict, shift=>result)
                next!(p)
            end
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
            exfiltrateAfter::Int=Inf,
            kws...)::AbstractDict

        safe_dict = ThreadSafeDict()
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
        if multi == :thread
            P = Progress(length(sets), desc=msg)
            Threads.@threads for (i, (shuffle,shift)) in sets
                data = shuffle_func(data, shuffle_pos...; shuffle_kws...)
                result = field.get_fields(σ(beh,shift), data; kws...)
                if postfunc != nothing
                    result = postfunc(result)
                end
                push!(safe_dict, shift=>result)
                next!(P)
                if mod(i, exfiltrateAfter)
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
                push!(safe_dict, (shift=shift,shuffle=shuffle)=>result)
            end
        elseif multi == :single
            @showprogress 0.1 msg for (i,(shuffle,shift)) in sets
                data = shuffle_func(data, shuffle_pos...; shuffle_kws...)
                result = field.get_fields(σ(beh,shift), data; kws...)
                if postfunc != nothing
                    result = postfunc(result)
                end
                push!(safe_dict, (shift=shift,shuffle=shuffle)=>result)
                if mod(i, exfiltrateAfter)
                    @exfiltrate
                end
            end
        else
            throw(ArgumentError("Unrecognized argument multi=$multi"))
        end
        safe_dict = Dict(safe_dict...)
        out = OrderedDict(key=>pop!(safe_dict, key) 
                          for key in sort([keys(safe_dict)...]))
        return out
    end


    function to_dataframe(shifts::AbstractDict; kws...) 
        table.to_dataframe(shifts, key_name="shift", kws...)
    end

    function info_to_dataframe(shifts::AbstractDict; kws...) where T <: AbstractArray
        table.to_dataframe(shifts, key_name="shift", name="info", kws...)
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


end
export timeshift
