__precompile__(false)
module Timeshift

    __revise_mode__ = :eval
    using DrWatson
    using Revise
    using Infiltrator
    using ProgressMeter, ProgressLogging
    using DataStructures: OrderedDict

    # Parent libary
    import ..GoalFetchAnalysis: Field
    import GoalFetchAnalysis.Field: adaptive, fixed
    import GoalFetchAnalysis.Field.preset: field_presets, return_preset_funcs
    import GoalFetchAnalysis.DI: Filt
    using DIutils
    Utils = DIutils # for compat with old checkpoint files

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
    export shift_func!, reset_shift!

    if !isdefined(Timeshift, :DictOfShiftOfUnit)
        DictOfShiftOfUnit{T<:Union{<:Real, NamedTuple}} =
                                                OrderedDict{T, OrderedDict} 
        AbsDictOfShiftOfUnit =
                         OrderedDict{<:Union{<:Real, NamedTuple}, OrderedDict} 
    end

    function isminutes(X::DataFrame)
        DIutils.dextrema(X.time)[1] < 1440.0 # assumes less than 24 hour recording
    end

    function ensureTimescale!(X::DataFrame; kws...)
        if isminutes(X; kws...)
            transform!(X, :time => (x->x.*60) => :time)
        else
            X
        end
    end

    
    function _functionalize(x::Union{Symbol,Function,Nothing})::Union{Nothing,Function}
        x !== nothing && x isa Symbol ? eval(x) : x
    end

    # -------------------- SHIFTING TYPES ---------------------------
    function shift_func!(data::DataFrame, shift::Real) 
        if :time_copy ∉ propertynames(data)
            data[!, :time_copy] = data[:, :time]
        end
        transform!(data, :time_copy => (t->t.+shift) =>:time)
    end
    const σ = shift_func!

    function reset_shift!(data::DataFrame) 
        if :time_copy ∈ propertynames(data)
            data[!, :time] = data[!, :time_copy]
        end
    end

    # -------------------- SHIFTED Receptive Fields --------------------------
    function get_field_shift(beh::DataFrame, data::DataFrame, shift::Real; 
            fieldfunc::Union{Function, Symbol}=Field.get_fields,
            postfunc::Union{Nothing,Function}=nothing,
            field_kws...)
        fieldfunc = fieldfunc isa Symbol ? eval(fieldfunc) : fieldfunc
        fieldobj  = fieldfunc(σ(beh, shift), data; field_kws...)
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
            result_dict::AbstractDict=OrderedDict(),
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
            grid=nothing,
            occ=nothing,
            fieldfunc::Union{Function, Symbol,Nothing}=nothing,
            gridfunc::Union{Function, Symbol,Nothing}=nothing,
            occfunc::Union{Function, Symbol,Nothing}=nothing,
            postfunc::Union{Function,Nothing}=nothing,
            metricfuncs::Union{Function,Symbol,
                         Vector{Symbol},Vector{Function},Nothing}=nothing,
            multi::Union{Bool, Symbol}=true,
            result_dict::AbstractDict=Timeshift.DictOfShiftOfUnit{Float64}(),
            progress::Bool=true,
            prog_fields::Bool=false,
            shiftbeh::Bool=true,
            thread_field::Bool=adaptive.thread_field_default,
            thread_fields::Bool=adaptive.thread_fields_default,
            overwrite_precache::Bool=false,
            grid_kws...)::OrderedDict

        @info "Timeshift" prog_fields

        # Process preset options
        if fieldpreset !== nothing
            f, g, o, m, p = return_preset_funcs(fieldpreset)
            fieldfunc    = fieldfunc   === nothing ? f : fieldfunc
            gridfunc     = gridfunc    === nothing ? g : gridfunc
            occfunc      = occfunc     === nothing ? o : occfunc
            metricsfuncs = metricfuncs === nothing ? m : metricfuncs
            postfunc     = postfunc    === nothing ? p : postfunc
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
        @info metricfuncs

        if filters !== nothing
            if Filt.filters_use_precache(filters) &&
                (overwrite_precache || Filt.missing_precache_input_cols(data,
                                                                 filters))
                @debug "precaching"
                Filt.precache!(data, beh, filters)
            end
            @debug "filtering"
            beh  = filtreg.filter(beh; filters, filter_skipmissingcols=true)[1]
        end
        grid = grid === nothing ? gridfunc(beh, props; grid_kws...) : grid
        occ  = occ  === nothing ? occfunc(beh, grid) : occ

        #ensureTimescale!(data)
        #ensureTimescale!(beh)
        shiftbeh ? reset_shift!(beh) : reset_shift!(data)

        if !(isdefined(Main, :PlutoRunner)) && progress
            prog = Progress(length(shifts), desc="Field shift calculations")
            prog.showspeed = true
        end
        @info shiftbeh
        if typeof(result_dict) <: Timeshift.AbsDictOfShiftOfUnit 
            result_dict = Timeshift.DictOfShiftOfUnit{keytype(result_dict)}(result_dict)
        end
        sizehint!(result_dict, length(shifts))
        for shift in shifts
            if shiftbeh
                @debug "shifting beh"
                shift = shift == 0 ? shift : shift * (-1)
                beh= σ(beh, shift) 
                @debug "register"
                beh, data = filtreg.register(beh, data; on="time",
                                              transfer=grid.props)
            else
                @debug "Shifting spike"
                data = σ(data, shift)
                @debug "register"
                beh, data = filtreg.register(beh, data; on="time",
                                              transfer=grid.props)
            end
            if shift ∈ keys(result_dict)
                continue
            end
            @debug "fieldfunc" thread_field thread_fields
            if postfunc === nothing
                push!(result_dict, shift =>
                            fieldfunc(data, grid, occ;
                                      splitby, metrics=metricfuncs,
                                      prog_fields,
                                      filters=nothing, # already filtered
                                      thread_field, thread_fields)
                           )
            else
                push!(result_dict, shift =>
                            postfunc(fieldfunc(data, grid, occ;
                                     splitby, metrics=metricfuncs,
                                     filters=nothing,
                                     thread_field, thread_fields))
                           )
            end
            #@infiltrate all(isnan.([value[:maxrate] for value in values(result_dict[shift])]))
            if !(isdefined(Main, :PlutoRunner)) && progress
                next!(prog)
            end
        end
        if !(isdefined(Main, :PlutoRunner)) && progress
            finish!(prog)
        end
        shiftbeh ? reset_shift!(beh) : reset_shift!(data)
        return result_dict
    end

    include("Timeshift/checkpoint.jl")
    include("Timeshift/dataframe.jl")
    include("Timeshift/operation.jl")
    include("Timeshift/shuffle.jl")
    include("Timeshift/crossval.jl")
    include("Timeshift/Keys.jl")
    include("Timeshift/types.jl")
    include("Timeshift/shiftmetrics.jl")

end
