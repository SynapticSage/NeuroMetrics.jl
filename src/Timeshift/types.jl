module types

    import Field
    import Table: to_dataframe
    import ..Timeshift: AbsDictOfShiftOfUnit, DictOfShiftOfUnit

    using DataFrames
    using DataStructures: OrderedDict

    export ShiftedField, ShiftedFields
    export getindex, getshifts, getunits

    #function Base.get(S::AbsDictOfShiftOfUnit)
    #    S.data
    #end
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
            k = collect(keys(s))
            kk = [k[1] for k in Tuple.(k)]
            i = pop!(index)
            select = kk .== i
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
        key = first(keys(S))
        collect(keys(S[key]))
    end

    mutable struct ShiftedField
        values::Vector{<:Field.ReceptiveField}
        keys::Vector{<:Real}
        shifts::Vector{<:Real}
        metrics::DataFrame

        function ShiftedField(data::OrderedDict)
            shifts = collect(keys(data))
            metrics = OrderedDict(k=>v.metrics for (k,v) in data)
            metrics = to_dataframe(metrics; key_name="shift")
            new(collect(values(data)), shifts, shifts, metrics)
        end

        function ShiftedField(data::OrderedDict, shifts::Vector{<:Real}) 
            metrics = OrderedDict(k=>v.metrics for (k,v) in data)
            metrics = to_dataframe(metrics; key_name="shift")
            new(collect(values(data)), shifts, shifts, metrics)
        end
    end

    mutable struct ShiftedFields
        values::Vector{ShiftedField}
        keys::Vector{NamedTuple}
        metrics::DataFrame

        function ShiftedFields(S::AbsDictOfShiftOfUnit)
            shifts = getshifts(S)
            units  = getunits(S)
            fields = OrderedDict{NamedTuple, ShiftedField}()
            metrics = Dict()
            for unit in units
                SF = ShiftedField(OrderedDict(shift=>get(S, shift, unit) for shift in shifts))
                fields[unit] = SF
                metrics[unit] = SF.metrics
            end
            new(collect(values(fields)), collect(keys(fields)), to_dataframe(metrics))
        end
    end

end
