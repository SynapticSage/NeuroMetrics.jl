module timeshift
    using DataFrames
    import ..field
    export get_field_shift
    export field
    using ThreadSafeDicts
    using DataStructures

    shift_func(data, shift) = transform(data, :time => (t->t.+shift) =>:time)
    const σ = shift_func

    function get_field_shift(beh::DataFrame, data::DataFrame, shift::Real; kws...)
        fieldobj = field.get_fields(σ(beh, shift), data; kws...)
    end

    function get_field_shift(beh::DataFrame, data::DataFrame,
            shifts::Union{StepRangeLen,Vector{Real}}; kws...)

        out = ThreadSafeDict()
        Threads.@threads for shift in shifts
            push!(out, shift=>field.get_fields(σ(beh,shift), data; kws...))
        end
        return out
    end

end
export timeshift
