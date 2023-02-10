module legacy

    using DataFrames
    using DataStructures: OrderedDict
    using ThreadSafeDicts
    using ProgressMeter

    import ..Timeshift: σ
    import ..Timeshift: Field

    """
        get_field_shift

    old way of getting shifted fields
    """
    function get_field_shift(beh::DataFrame, data::DataFrame,
            shifts::Union{StepRangeLen,Vector{T}} where T <: Real; 
            fieldfunc::Union{Function, Symbol}=Field.legacy.get_fields,
            postfunc::Union{Function,Nothing}=nothing,
            as::Type=OrderedDict,
            multi::Union{Bool, Symbol}=true,
            safe_dict::AbstractDict=ThreadSafeDict(),
            squeeze_unity::Bool=false,
            kws...)

        if multi isa Bool
            multi = multi ? :thread : :single
        end
        @info "Starting multi=$multi"

        p = Progress(length(shifts), desc="Field shift calculations")
        if multi == :thread
            Threads.@threads for shift in shifts
                if shift ∈ keys(safe_dict)
                    continue
                end
                result = fieldfunc(σ(beh, shift), data; kws...)
                if postfunc !== nothing
                    result = postfunc(result)
                end
                push!(safe_dict, shift=>result)
                next!(p)
            end
        elseif multi == :single
            @showprogress 0.1 "shifts=$shifts" for shift in shifts
                if shift ∈ keys(safe_dict)
                    continue
                end
                result = fieldfunc(σ(beh,shift), data; kws...)
                if postfunc !== nothing
                    result = postfunc(result)
                end
                push!(safe_dict, shift=>result)
                next!(p)
            end
        end
        safe_dict = Dict(safe_dict...)
        out = as(key=>pop!(safe_dict, key) for key in sort([keys(safe_dict)...]))
        if length(keys(safe_dict)) == 1 && squeeze_unity
            out = out[first(keys(out))]
        end
        return out
    end

    # -------------------- SHIFTED Receptive Fields --------------------------
    function get_field_shift(beh::DataFrame, data::DataFrame, shift::Real; 
            fieldfunc::Union{Function, Symbol}=Field.legacy.get_fields,
            postfunc::Union{Nothing,Function}=nothing,
            field_kws...)
        fieldfunc = fieldfunc isa Symbol ? eval(fieldfunc) : fieldfunc
        fieldobj = fieldfunc(σ(beh, shift), data; field_kws...)
        postfunc !== nothing ? fieldobj = postfunc(fieldobj) : fieldobj
    end

end
