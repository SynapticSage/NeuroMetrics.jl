module Timeshift

    __revise_mode__ = :eval
    using DrWatson
    using Revise
    using Infiltrator
    using ProgressMeter
    using DataStructures: OrderedDict

    # Parent libary
    import Field
    import Load: register
    adaptive = Field.adaptive

    # Julia packages
    using DataStructures
    using ThreadSafeDicts
    using DataFrames, DataFramesMeta
    using Distributions
    using NaNStatistics

    macro identity(n)
        return n
    end

    # Exports
    export Field
    export get_field_shift
    export shifted_fields

    # -------------------- SHIFTING TYPES ---------------------------
    shift_func(data::DataFrame, shift::Real) = 
             transform(data, :time => (t->t.+shift) =>:time, copycols=false)
    const σ = shift_func

    # -------------------- SHIFTED Receptive Fields --------------------------
    function get_field_shift(beh::DataFrame, data::DataFrame, shift::Real; 
            fieldfunc::Union{Function, Symbol}=Field.get_fields,
            postfunc::Union{Nothing,Function}=nothing,
            field_kws...)
        fieldfunc = fieldfunc isa Symbol ? eval(fieldfunc) : fieldfunc
        fieldobj = fieldfunc(σ(beh, shift), data; field_kws...)
        postfunc !== nothing ? fieldobj = postfunc(fieldobj) : fieldobj
    end

    """
        get_field_shift

    old way of getting shifted fields
    """
    function get_field_shift(beh::DataFrame, data::DataFrame,
            shifts::Union{StepRangeLen,Vector{T}} where T <: Real; 
            fieldfunc::Union{Function, Symbol}=Field.get_fields,
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

    """
        function shifted_fields(beh::DataFrame, data::DataFrame,
                shifts::Union{StepRangeLen,Vector{T}} where T <: Real,
                props::Vector; 
                splitby::Vector=[:unit],
                fieldfunc::Union{Function, Symbol}=adaptive.ulanovsky,
                gridfunc::Union{Function, Symbol}=adaptive.get_grid,
                occfunc::Union{Function, Symbol}=adaptive.get_occupancy,
                postfunc::Union{Function,Nothing}=nothing,
                multi::Union{Bool, Symbol}=true,
                safe_dict::AbstractDict=ThreadSafeDict(),
                grid_kws...)::OrderedDict

    Newer versioin of get_field_shift that works with new adaptive and fixed
    types
    """
    function shifted_fields(beh::DataFrame, data::DataFrame,
            shifts::Union{StepRangeLen,Vector{T}} where T <: Real,
            props::Vector; 
            splitby::Vector=[:unit],
            fieldfunc::Union{Function, Symbol}=adaptive.ulanovsky,
            gridfunc::Union{Function, Symbol}=adaptive.get_grid,
            occfunc::Union{Function, Symbol}=adaptive.get_occupancy,
            postfunc::Union{Function,Nothing}=nothing,
            metrics::Union{Function,Symbol,
                           Vector{Symbol},Vector{Function},Nothing}=nothing,
            multi::Union{Bool, Symbol}=true,
            result_dict::AbstractDict=ThreadSafeDict(),
            grid_kws...)::OrderedDict

        if metrics isa Symbol
            metrics = [eval(metrics)]
        elseif metrics isa Function
            metrics = [metrics]
        elseif metrics isa Vector{Symbol}
            metrics = [eval(metric) for metric in metrics]
        end

        if multi isa Bool
            multi = multi ? :thread : :single
        end

        grid = gridfunc(beh, props; grid_kws...)
        occ  = occfunc(beh, grid)
        data = dropmissing(data, grid.props)

        if maximum(abs.(shifts)) >= 0.25 && median(diff(beh.time)) < (1/30)
            shift_correct = true
            @info "correcting shifts for minutes"
            shifts = collect(shifts)
            shifts ./= 60
        else
            shift_correct = false
        end

        prog = Progress(length(shifts), desc="Field shift calculations")
        for shift in shifts
            if shift ∈ keys(result_dict)
                continue
            end
            _, data = register(σ(beh, shift), data; on="time",
                               transfer=grid.props)
            result = fieldfunc(data, grid, occ; splitby)
            if postfunc !== nothing
                result = postfunc(result)
            end
            shift = shift_correct ? shift*60 : shift
            push!(result_dict, shift=>result)
            if metrics !== nothing
                for (unit,metric) in Iterators.product(keys(result), metrics)
                    push!(result[unit].metrics, Symbol(metric)=>metric(result[unit]))
                end
            end
            next!(prog)
        end
        result_dict = Dict(result_dict...)
        @infiltrate
        out = DictOfShiftOfUnit(
                          key=>pop!(result_dict, key) 
                          for key in sort([keys(result_dict)...])
                         )
        return out
    end

    DictOfShiftOfUnit = OrderedDict{Union{<:Real,<:NamedTuple}, Dict}
    function Base.get(S::DictOfShiftOfUnit)
        S.data
    end
    function Base.get(S::DictOfShiftOfUnit, shift::Float64, index...) 
        s = S[shift]
        while !(isempty(index))
            k = collect(keys(s))
            kₙ = [k[1] for k in Tuple.(k)]
            select = kₙ .== index[1]
            if !(any(select))
                @error "No matches"
            end
            s = s[k[select][1]]
        end
        return s
    end
    function Base.get(S::DictOfShiftOfUnit, shift::Float64,
             index::NamedTuple) 
    end

    mutable struct ShiftedField
        data::OrderedDict{<:Real,<:Field.ReceptiveField}
        shifts::Vector{<:Real}
        function ShiftedField(data::DictOfShiftOfUnit)
            shifts = keys(data)
            new(data, shifts)
        end
        ShiftedField(data::OrderedDict, shifts::Vector{<:Real}) = new(data, 
                                                                      shifts)
    end

    mutable struct ShiftedFields
        data::AbstractDict{<:NamedTuple, <:ShiftedField}
        function ShiftedFields(S::DictOfShiftOfUnit)
            vals = values(S.data)
            vals = [ShiftedField(val) for val in vals]
            kys  = keys(S.data)
            new(OrderedDict(ky=>val for (ky,val) in zip(kys,vals)))
        end
    end

    #@recipte plot_shiftedfields(S::ShiftedFields, shifts::Int, neurons::Int)
    #end

    using Reexport
    include(srcdir("Timeshift", "checkpoint.jl"))
    @reexport using .checkpoint
    include(srcdir("Timeshift", "dataframe.jl"))
    @reexport using .dataframe
    include(srcdir("Timeshift", "operation.jl"))
    @reexport using .operation
    include(srcdir("Timeshift", "shuffle.jl"))
    @reexport using .shuffle
    include(srcdir("Timeshift", "crossval.jl"))
    @reexport using .crossval
    include(srcdir("Timeshift", "Keys.jl"))
    @reexport using .Keys

end
