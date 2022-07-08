module Timeshift

    __revise_mode__ = :eval
    using DrWatson
    using Revise
    using Infiltrator
    using ProgressMeter
    using DataStructures: OrderedDict

    # Parent libary
    import Field
    import Load: register, filterAndRegister
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
    export getshifts, getunits
    export ShiftedField, ShiftedFields, DictOfShiftOfUnit

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
            filters::Union{OrderedDict, Nothing}=nothing,
            splitby::Vector=[:unit],
            fieldfunc::Union{Function, Symbol}=adaptive.ulanovsky,
            gridfunc::Union{Function, Symbol}=adaptive.get_grid,
            occfunc::Union{Function, Symbol}=adaptive.get_occupancy,
            postfunc::Union{Function,Nothing}=nothing,
            metricfuncs::Union{Function,Symbol,
                           Vector{Symbol},Vector{Function},Nothing}=nothing,
            multi::Union{Bool, Symbol}=true,
            result_dict::AbstractDict=ThreadSafeDict(),
            grid_kws...)::OrderedDict

        if metricfuncs isa Symbol
            metricfuncs = [eval(metricfuncs)]
        elseif metrics isa Function
            metricfuncs = [metricfuncs]
        elseif metrics isa Vector{Symbol}
            metricfuncs = [eval(metricfuncs) for metricfuncs in metricfuncs]
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
            _, data = filterAndRegister(σ(beh, shift), data; on="time",
                               transfer=grid.props, filters)
            result = fieldfunc(data, grid, occ; splitby)
            if postfunc !== nothing
                result = postfunc(result)
            end
            shift = shift_correct ? shift*60 : shift
            push!(result_dict, shift=>result)
            if metricfuncs !== nothing
                for (unit,metricfuncs) in Iterators.product(keys(result), metricfuncs)
                    name, calculation = Symbol(metricfuncs), metricfuncs(result[unit])
                    push!(result[unit].metricfuncs, name=>calculation)
                end
            end
            next!(prog)
        end
        result_dict = Dict(result_dict...)
        out = OrderedDict(
                          key=>pop!(result_dict, key) 
                          for key in sort([keys(result_dict)...])
                         )
        out = Timeshift.DictOfShiftOfUnit{keytype(out)}(out)
        return out
    end

    DictOfShiftOfUnit{T<:Union{<:Real, NamedTuple}} =
                                            OrderedDict{T, OrderedDict} 
    AbsDictOfShiftOfUnit =
                     OrderedDict{<:Union{<:Real, NamedTuple}, OrderedDict} 

    #function Base.get(S::AbsDictOfShiftOfUnit)
    #    S.data
    #end
    function Base.get(S::T where T<:AbsDictOfShiftOfUnit, shift::Float64, index...) 
        s = S[shift]
        index = [index...]
        while !(isempty(index))
            k = collect(keys(s))
            kk = [k[1] for k in Tuple.(k)]
            select = kk .== pop!(index)
            if !(any(select))
                @error "No matches"
            end
            s = s[k[select][1]]
        end
        return s
    end
    Base.get(S::AbsDictOfShiftOfUnit, shift::Float64, index::NamedTuple) = 
                                                            S[shift][index]
    function getshifts(S::AbsDictOfShiftOfUnit)
        collect(keys(S))
    end
    function getunits(S::AbsDictOfShiftOfUnit)
        collect(keys(S[1]))
    end

    mutable struct ShiftedField
        data::OrderedDict{<:Real,<:Field.ReceptiveField}
        shifts::Vector{<:Real}
        metrics::DataFrame
        function ShiftedField(data::OrderedDict)
            shifts = collect(keys(data))
            new(data, shifts)
        end
        ShiftedField(data::OrderedDict, shifts::Vector{<:Real}) = new(data, 
                                                                      shifts)
    end

    mutable struct ShiftedFields
        data::AbstractDict{<:NamedTuple, <:ShiftedField}
        metrics::DataFrame
        function ShiftedFields(S::AbsDictOfShiftOfUnit)
            shifts = getshifts(S)
            units  = getunits(S)
            fields = OrderedDict{NamedTuple, ShiftedField}()
            for unit in units
                @infiltrate
                SF = ShiftedField(OrderedDict(shift=>get(S, shift, unit) for shift in shifts))
                fields[unit] = SF
            end
            new(fields)
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
