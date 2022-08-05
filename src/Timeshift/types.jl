module types

    import Field
    import Field.metrics: metric_ban, apply_metric_ban, unstackMetricDF
    import Table: to_dataframe, vec_arrayofarrays!
    import ..Timeshift: AbsDictOfShiftOfUnit, DictOfShiftOfUnit
    import Utils

    using DataFrames
    using DataStructures: OrderedDict
    using Missings
    using AxisArrays
    using Infiltrator

    export ShiftedField, ShiftedFields
    export getshifts, getunits

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
    # Cache version
    function cachespikecount(spikes::DataFrame)
          groupby_summary_condition_column(spikes, :unit,
                                  x->x.spikecount .> req_spikes, # > 50 spikes
                                  nrow=>:spikecount; name=:spikecount_cache)
    end
    spikecountcached = OrderedDict(:spikecount_cache => x -> x) # stored true or false

    # Cache version
    function cachetrajdiversity(spikes::DataFrame)
        groupby_summary_condition_column(spikes, :unit,
                                x -> x.percount.>=req_traj, 
                                :period => (x->length(unique(x))) =>
                                :trajdiversity; name=:trajdiversity_cache)
    end
    trajdiversitycached       = OrderedDict(:trajdiversity_cache =>
                                         x -> x) # stored true or false

    """
    stores functions used to precache each field
    """
    precache_funcs = Dict(:trajdiversity_cache => cachetrajdiversity, 
                         :spikecount_cache     => cachespikecount)
    """
        precache_field_reqs

    stores the fields required to precache
    """
    precache_field_reqs = Dict(:trajdiversity_cache => [:period])

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

        function ShiftedField(data::OrderedDict{<:Any, <:Field.ReceptiveField})
            shifts = collect(keys(data))
            metrics = OrderedDict(k=>apply_metric_ban(v.metrics) 
                                  for (k,v) in data 
                                  if !(ismissing(v)))
            metricDF = to_dataframe(metrics; key_name="shift", explode=false)
            metricDF = unstackMetricDF(metricDF)
            new(collect(values(data)), shifts, shifts, metricDF)
        end

        function ShiftedField(data::OrderedDict, shifts::Vector{<:Real}) 
            metrics = OrderedDict(k=>apply_metric_ban(v.metrics) for (k,v) in data
                                  if !(ismissing(v)))
            metricsDF = to_dataframe(metrics; key_name="shift", explode=false)
            metricDF  = unstackMetricDF(metricDF)
            new(collect(values(data)), shifts, shifts, metrics)
        end
    end
    OrderedDict(SF::ShiftedField) = OrderedDict(zip(SF.keys, SF.values))
    Base.Dict(SF::ShiftedField)   = OrderedDict(zip(SF.keys, SF.values))
    to_dataframe(SF::ShiftedField; kws...) = to_dataframe(Dict(SF); 
                                                          key_name=["shift","property"],
                                                         kws...)
    function Base.getindex(SF::ShiftedField, shift::T where T<:Real)
        SF.values[findfirst(SF.keys .== shift)]
    end

    

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

    function getunits(sfs::ShiftedFields)
        sfs.keys
    end
    function getshifts(sfs::ShiftedFields)
        sfs[first(sfs.keys)].keys
    end

    function Base.getindex(SFs::ShiftedFields, unit::Real) 
        D = OrderedDict(SFs)
        D[keymatch_topcomponent(D, unit)]
    end
    function Base.getindex(SFs::ShiftedFields, unit::NamedTuple) 
        SFs.values[findfirst([tuple(k) == tuple(unit) for k in SFs.keys])]
    end
    function Base.getindex(SFs::ShiftedFields, unit::NamedTuple, shift::Real) 
        SFs[unit][shift]
    end
    Base.getindex(SFs::ShiftedFields, unit::Real, shift::Real) = SFs[unit][shift]

    OrderedDict(SF::ShiftedFields)   = OrderedDict(zip(SF.keys, SF.values))
    Base.Dict(SF::ShiftedFields)     = OrderedDict(zip(SF.keys, SF.values))
    to_dataframe(SFs::ShiftedFields; kws...) = to_dataframe(Dict(SFs);
                                                           kws...)

    """
    Creates a matrix of shifted field objects
    """
    function matrixform(sfs::ShiftedFields)::AxisArray
        units, shifts = getunits(sfs), getshifts(sfs)
        M = Array{Field.ReceptiveField, 2}(undef, length(units), length(shifts))
        for (u, unit) in enumerate(units)
            for (s, shift) in enumerate(shifts)
                val = sfs[unit, shift]
                M[u,s] = val
            end
        end
        units = [unit[1] for unit in units]
        AxisArray(M, Axis{:unit}(units), Axis{:shift}(shifts))
    end

    """
    Creates a tensor of shifted field data
    """
    function tensorform(sfs::ShiftedFields)::AxisArray
        M = matrixform(sfs)
        grid = M[1,1].grid
        val = M[1,1].rate
        selector = [Utils.na, Utils.na, (Colon() for i in 1:ndims(val))...]
        results = []
        for mm in eachrow(M)
            res = cat([getproperty(m,:rate) for m in mm]...; dims=4)
            push!(results,res)
        end
        results = cat(results...; dims=3)
        results = permutedims(results, (3,4,1,2))
        fieldaxes = [Axis{Symbol(grid.props[i])}(grid.centers[i]) for
                     i in 1:length(grid.props)]
        AxisArray(results, M.axes..., fieldaxes...)
    end

end
