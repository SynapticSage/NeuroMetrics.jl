module convert

    export to_dataframe

    """
    to_dataframe

    #purpose
    converts a field dictionary (multiple units) into a field dataframe
    """
    #function to_dataframe(fields::AbstractDict; kws...)
    #    fields = Dict(key=>fields[key] for key in keys(fields))
    #end
    function to_dataframe(fields::AbstractDict; other_labels=Dict(), 
            key_name::Union{Nothing,String}=nothing, kws...)
        D = DataFrame()
        for (key, field) in fields
            @debug "key=$key"
            if key_name != nothing #&& (key isa NamedTuple || key isa Dict)
                other_labels[key_name] = key
            elseif (key isa NamedTuple) || (key isa Dict)
                key_dict = Dict(string(k)=>v for (k, v) in pairs(key))
                other_labels = merge(other_labels, key_dict)
            else
                @warn "unhandled key_name=$key"
                other_labels["unnamed"] = key
            end
            if fields[key] != nothing
                append!(D, to_dataframe(fields[key];
                                        other_labels=other_labels, kws...))
            end
        end
        return D
    end
    function to_dataframe(fields::T where T <: NamedTuple; kws...)
        D = to_dataframe(Utils.namedtuple_to_dict(fields); kws...) 
        return D
    end

    """
    field_to_dataframe(field::AbstractArray; other_labels=Dict())

    +purpose: converts a single field matrix into a field dataframe
    """
    function to_dataframe(F::Union{AbstractArray,Real};
            props::Vector{String}=Vector{String}([]), grid::Tuple=(),
            gridtype="center", other_labels=Dict(), name::String="value",
            explode::Bool=true)
        D = ndgrid((1:size(F,i) for i in 1:ndims(F))...)
        D = OrderedDict{String,Any}("dim_$d"=>vec(D[d])
                              for d in 1:ndims(F))
        if typeof(F) <: Real
            F = [F]
        end
        if ~isempty(props)
            if gridtype == "edge"
                grid = edge_to_center.(grid)
            end
            grid = ndgrid(grid...)
        end
        _clean_label_values(other_labels)
        for (label, value) in other_labels
            if label isa Symbol
                label = String(label)
            end
            D[label] = value
        end
        if explode
            for (prop, G) in zip(props, grid)
                D[prop] = vec(G)
            end
            D[name] = vec(F);
        else
            for (key, value) in D
                D[key] = [value]
            end
            D[name] = [F,];
        end
        D = DataFrame(D)
    end

    function _clean_label_values(L::AbstractDict)
        for (k,v) in L
            typev = typeof(v)
            if typev <: AbstractVector
                if eltype(v) == String
                    v = join(v, "-")
                else
                    v = string(v)
                end
                L[k] = v
            end
        end
    end
end
