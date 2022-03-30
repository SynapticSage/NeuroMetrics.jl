module timeshift
    using ..field
    using Base.Threads.@spawn
    function get_field_shift(beh::DataFrame, data::DataFrame, shift::Real; kws...)
        beh = copy(beh)
        beh.time = beh.time + shift
        fieldobj = get_fields(beh, data)
    end
    function get_field_shift(beh::DataFrame, data::DataFrame, shifts::Vector{Real},
            kws...)
        (@spawn get_fields(beh, data; kws...) for shift in shifts)
    end
end
