module Timeshift

    __revise_mode__ = :eval
    using DrWatson
    using Revise
    using Infiltrator
    using ProgressMeter
    using DataStructures: OrderedDict

    # Parent libary
    import Field
    import Field: adaptive, fixed
    import Load: register, filterAndRegister
    import Field.preset: field_presets, return_preset_funcs
    import Table: to_dataframe
    import Munge.chrono: ensureTimescale

    # Julia packages
    using DataStructures
    using ThreadSafeDicts
    using DataFrames, DataFramesMeta
    using Distributions
    using NaNStatistics


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
            shifts::Union{StepRangeLen,StepRange,Vector{T}} where T <: Real,
            props::Vector; 
            filters::Union{OrderedDict, Nothing}=nothing,
            splitby::Vector=[:unit],
            fieldpreset::Union{Nothing,Symbol}=:yartsev,
            fieldfunc::Union{Function, Symbol,Nothing}=nothing,
            gridfunc::Union{Function, Symbol,Nothing}=nothing,
            occfunc::Union{Function, Symbol,Nothing}=nothing,
            postfunc::Union{Function,Nothing}=nothing,
            metricfuncs::Union{Function,Symbol,
                           Vector{Symbol},Vector{Function},Nothing}=nothing,
            multi::Union{Bool, Symbol}=true,
            result_dict::AbstractDict=ThreadSafeDict(),
            progress::Bool=true,
            grid_kws...)::OrderedDict

        # Process preset options
        if fieldpreset !== nothing
            fieldfunc, gridfunc, occfunc, metricfuncs, postfunc =
                return_preset_funcs(fieldpreset)
        else
            fieldfunc   = fieldfunc isa Symbol ? eval(fieldfunc) : fieldfunc
            occfunc     = occfunc isa Symbol ? eval(occfunc) : occfunc
            gridfunc    = gridfunc isa Symbol ? eval(gridfunc) : gridfunc
            postfunc    = postfunc isa Symbol ? eval(postfunc) : postfunc
        end
        if metricfuncs isa Symbol
            metricfuncs = [eval(metricfuncs)]
        elseif metricfuncs isa Function
            metricfuncs = [metricfuncs]
        elseif metricfuncs isa Vector{Symbol}
            metricfuncs = [eval(metricfuncs) for metricfuncs in metricfuncs]
        end
        if multi isa Bool
            multi = multi ? :thread : :single
        end

        grid = gridfunc(beh, props; grid_kws...)
        occ  = occfunc(beh, grid)
        data = dropmissing(data, grid.props)

        data = ensureTimescale(data)
        beh  = ensureTimescale(beh)

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
            push!(result_dict, shift=>result)
            if metricfuncs !== nothing
                for (unit,metricfuncs) in Iterators.product(keys(result), metricfuncs)
                    name, calculation = Symbol(metricfuncs), metricfuncs(result[unit])
                    push!(result[unit].metrics, name=>calculation)
                end
            end
            if !(isdefined(Main, :PlutoRunner)) && progress
                next!(prog)
            end
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
    function Base.get(S::T where T<:AbsDictOfShiftOfUnit, shift::T where T<:Real, index...) 
        s = S[shift]
        getindex(s, index...)
    end
    function getindex(s::AbstractDict{NamedTuple,<:Any}, index...)
        index = [index...]
        while !(isempty(index))
            k = collect(keys(s))
            kk = [k[1] for k in Tuple.(k)]
            select = kk .== pop!(index)
            #@info "select = $select"
            if !(any(select))
                @error "No matches"
            end
            s = s[k[select][1]]
        end
        return s
    end
    function Base.get(S::T where T<:AbsDictOfShiftOfUnit, shift::Colon, index...) 
        OrderedDict(k=>getindex(v,index...) for (k,v) in S)
    end
    function Base.get(S::T where T<:AbsDictOfShiftOfUnit, shift::T where T<:Real) 
        S[shift]
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
            metrics = OrderedDict(k=>v.metrics for (k,v) in data)
            metrics = to_dataframe(metrics)
            new(data, shifts, metrics)
        end
        function ShiftedField(data::OrderedDict, shifts::Vector{<:Real}) 
            metrics = OrderedDict(k=>v.metrics for (k,v) in data)
            metrics = to_dataframe(metrics)
            new(data, shifts, metrics)
        end
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
