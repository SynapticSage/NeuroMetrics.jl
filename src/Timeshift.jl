module Timeshift

    __revise_mode__ = :eval
    using DrWatson
    using Revise
    using Infiltrator
    using ProgressMeter, ProgressLogging
    using DataStructures: OrderedDict

    # Parent libary
    import Field
    import Field: adaptive, fixed
    import Load: utils
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


    
    function _functionalize(x::Union{Symbol,Function,Nothing})::Union{Nothing,Function}
        x !== nothing && x isa Symbol ? eval(x) : x
    end

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
        shifted_fields(beh::DataFrame, data::DataFrame,
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

    Newer versioin of `get_field_shift` that works with new adaptive and fixed
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
        end
        fieldfunc   = _functionalize(fieldfunc)
        occfunc     = _functionalize(occfunc)
        gridfunc    = _functionalize(gridfunc)
        postfunc    = _functionalize(postfunc)
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

        beh  = filters !== nothing ? utils.filter(beh; filters, filter_skipmissingcols=true)[1] : beh
        grid = gridfunc(beh, props; grid_kws...)
        occ  = occfunc(beh, grid)

        data = ensureTimescale(data)
        beh  = ensureTimescale(beh)

        if !(isdefined(Main, :PlutoRunner)) && progress
            prog = Progress(length(shifts), desc="Field shift calculations")
        end
        for shift in shifts
            if shift ∈ keys(result_dict)
                continue
            end
            shift_beh    = σ(beh, shift)
            _, data_filt = utils.register(shift_beh, data; on="time",
                               transfer=grid.props)
            data_filt = dropmissing(data_filt, grid.props)
            result    = fieldfunc(data_filt, grid, occ; splitby)
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
    include(srcdir("Timeshift", "types.jl"))
    @reexport using .types

end
