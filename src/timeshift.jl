module timeshift
    using DataFrames
    import ..field
    export get_field_shift

    function get_field_shift(beh::DataFrame, data::DataFrame, shift::Real; kws...)
        beh = copy(beh)
        beh.time = beh.time .+ shift
        fieldobj = field.get_fields(beh, data; kws...)
    end

    function get_field_shift(beh::DataFrame, data::DataFrame,
            shifts::Vector{Real}; kws...)

        out = []
        Threads.@threads for shift in shifts
            append!(out, get_fields(beh, data; kws...) for shift in shifts)
        end
    end

end
export timeshift
