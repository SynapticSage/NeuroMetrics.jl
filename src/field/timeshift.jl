module timeshift
    using DataFrames
    import ..field
    export get_field_shift
    export field
    using ThreadSafeDicts
    using DataStructures
    include("table.jl")
    export table

    # -------------------- SHIFTING TYPES ---------------------------
    shift_func(data, shift) = transform(data, :time => (t->t.+shift) =>:time)
    const σ = shift_func

    # -------------------- SHIFTED Receptive Fields --------------------------
    function get_field_shift(beh::DataFrame, data::DataFrame, shift::Real; kws...)
        fieldobj = field.get_fields(σ(beh, shift), data; kws...)
    end

    function get_field_shift(beh::DataFrame, data::DataFrame,
            shifts::Union{StepRangeLen,Vector{Real}}; 
            multithread::Bool=true,
            postfunc::Union{Function,Nothing}=nothing,
            kws...)

        safe_dict = ThreadSafeDict()
        if multithread
            Threads.@threads for shift in shifts
                result = field.get_fields(σ(beh,shift), data; kws...)
                if postfunc
                    result = postfunc(result)
                end
                push!(safe_dict, shift=>result)
            end
        else
            for shift in shifts
                result = field.get_fields(σ(beh,shift), data; kws...)
                if postfunc
                    result = postfunc(result)
                end
                push!(safe_dict, shift=>result)
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

end
export timeshift
