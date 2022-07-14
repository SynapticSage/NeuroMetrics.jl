module types

    import Field
    import Table: to_dataframe
    import ..Timeshift: AbsDictOfShiftOfUnit, DictOfShiftOfUnit

    using DataFrames
    using DataStructures: OrderedDict
    using Missings
    using Infiltrator

    export ShiftedField, ShiftedFields
    export getindex, getshifts, getunits

    function Base.get(S::T where T<:AbsDictOfShiftOfUnit, shift::T where T<:Real, index...) 
        s = S[shift]
        getindex(s, index...)
    end
    function getindex(s::AbstractDict{NamedTuple,<:Any}, index::Union{NamedTuple}...)
        index = [index...]
        while !(isempty(index))
            i = pop!(index)
            if i ∉ keys(s)
                @error "No matches key=$i"
            end
            s = s[i]
        end
        return s
    end
    function getindex(s::AbstractDict{NamedTuple,<:Any}, index::Union{Tuple,Array,Real}...)
        index = [index...]
        while !(isempty(index))
            select =  keymatch_topcomponent(s, pop!(index))
            s = s[select]
        end
        return s
    end
    function keymatch_topcomponent(s::AbstractDict{NamedTuple, <:Any}, 
            i::Union{Tuple, Array, Real})
        k = collect(keys(s))
        kk = [k[1] for k in Tuple.(k)]
        select = kk .== i
        if !(any(select))
            @error "No matches"
        end
        k[select][1]
    end
    function Base.get(S::T where T<:AbsDictOfShiftOfUnit, shift::Colon, index...) 
        OrderedDict(k=>getindex(v,index...) for (k,v) in S)
    end
    function Base.get(S::T where T<:AbsDictOfShiftOfUnit, shift::T where T<:Real) 
        S[shift]
    end
    function Base.get(S::AbsDictOfShiftOfUnit, shift::Float64, index::NamedTuple)
        if index ∈ keys(S[shift])
            S[shift][index]
        else
            missing
        end
    end
    function getshifts(S::AbsDictOfShiftOfUnit)
        collect(keys(S))
    end
    function getunits(S::AbsDictOfShiftOfUnit)
        key = first(keys(S))
        collect(keys(S[key]))
    end

    mutable struct ShiftedField
        values::Vector{Union{<:Field.ReceptiveField, Missing}}
        keys::Vector{<:Real}
        shifts::Vector{<:Real}
        metrics::DataFrame

        function ShiftedField(data::OrderedDict)
            shifts = collect(keys(data))
            metrics = OrderedDict(k=>v.metrics for (k,v) in data if !(ismissing(v)))
            metrics = to_dataframe(metrics; key_name="shift")
            new(collect(values(data)), shifts, shifts, metrics)
        end

        function ShiftedField(data::OrderedDict, shifts::Vector{<:Real}) 
            metrics = OrderedDict(k=>v.metrics for (k,v) in data if !(ismissing(v)))
            metrics = to_dataframe(metrics; key_name="shift")
            new(collect(values(data)), shifts, shifts, metrics)
        end
    end
    OrderedDict(SF::ShiftedField) = OrderedDict(zip(SF.keys, SF.values))
    Base.Dict(SF::ShiftedField)   = OrderedDict(zip(SF.keys, SF.values))
    to_dataframe(SF::ShiftedField; kws...) = to_dataframe(Dict(SF); 
                                                          key_name=["shift","property"],
                                                         kws...)

    mutable struct ShiftedFields
        values::Vector{Union{Missing,ShiftedField}}
        keys::Vector{NamedTuple}
        metrics::DataFrame

        function ShiftedFields(S::AbsDictOfShiftOfUnit)
            shifts = getshifts(S)
            units  = getunits(S)
            fields = OrderedDict{NamedTuple, Union{Missing,ShiftedField}}()
            metrics = Dict()
            for unit in units
                SF = OrderedDict(shift=>get(S, shift, unit) for shift in shifts)
                SF = ShiftedField(SF)
                fields[unit]  = SF
                metrics[unit] = SF.metrics
            end
            new(collect(values(fields)), collect(keys(fields)), to_dataframe(metrics))
        end
    end
    function Base.getindex(SFs::ShiftedFields, i::Real) 
        D = OrderedDict(SFs)
        D[keymatch_topcomponent(D, i)]
    end
    OrderedDict(SF::ShiftedFields)   = OrderedDict(zip(SF.keys, SF.values))
    Base.Dict(SF::ShiftedFields)     = OrderedDict(zip(SF.keys, SF.values))
    to_dataframe(SFs::ShiftedFields; kws...) = to_dataframe(Dict(SFs);
                                                           kws...)

end
