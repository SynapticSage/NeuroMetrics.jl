module timeshift
    using DataFrames
    import ..field
    export get_field_shift
    export field
    using ThreadSafeDicts
    using DataStructures
    using DataFramesMeta
    using ProgressMeter
    include("../table.jl")
    include("../shuffle.jl")
    export table

    # -------------------- SHIFTING TYPES ---------------------------
    shift_func(data, shift) = transform(data, :time => (t->t.+shift) =>:time)
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
            kws...)

        safe_dict = ThreadSafeDict()
        p = Progress(length(shifts), desc="field shift calculations")
        if multithread
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
            multithread::Bool=true,
            nShuffle::Int=0, shuffle_pos::Tuple=(), shuffle_kws::NamedTuple=(;), shuffle_func=shuffle.by,
            postfunc::Union{Function,Nothing}=nothing,
            kws...)

        safe_dict = ThreadSafeDict()
        sets = collect(product(1:nShuffle, shifts))
        p = Progress(length(sets), desc="SHUFFLED field shift calculations")
        if multithread
            Threads.@threads for (shuffle,shift) in sets
                shuffle_kws[:shuffledist_df] = beh
                spikes = shuffle_func(spikes; shuffle_kws...)
                result = field.get_fields(σ(beh,shift), data; kws...)
                if postfunc != nothing
                    result = postfunc(result)
                end
                push!(safe_dict, shift=>result)
            end
        else
            @showprogress for (shuffle,shift) in sets
                shuffle_kws[:shuffledist_df] = beh
                spikes = shuffle_func(spikes; shuffle_kws...)
                result = field.get_fields(σ(beh,shift), data; kws...)
                if postfunc != nothing
                    result = postfunc(result)
                end
                push!(safe_dict, (shift=shift,shuffle=shuffle)=>result)
            end
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
