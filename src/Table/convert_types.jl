module convert_types

    using DataFrames
    export to_dataframe
    using LazyGrids: ndgrid
    using DataStructures
    import Utils
    using Infiltrator

    """
    to_dataframe

    #purpose
    converts a field dictionary (multiple units) into a field dataframe
    """
    #function to_dataframe(fields::AbstractDict; kws...)
    #    fields = Dict(key=>fields[key] for key in keys(fields))
    #end
    function to_dataframe(fields::AbstractDict; other_labels=Dict(), 
            key_name::Union{Nothing,Symbol,String,Vector}=nothing,
            level::Int=0, level_names::Dict=Dict(), kws...)::DataFrame

        level += 1
        @debug level

        D = DataFrame()
        for key in keys(fields)

            if (key isa NamedTuple) || (key isa Dict)
                key_dict = Dict(string(k)=>v for (k, v) in pairs(key))
                other_labels = merge(other_labels, key_dict)

            elseif key_name !== nothing #&& (key isa NamedTuple || key isa Dict)

                if (typeof(key_name) <: Vector)
                    kn = pop!(key_name)
                elseif (key_name == :keyboard)
                    if level ∉ keys(level_names)
                        println("Name key for key=$key, typeof(key)=$key")
                        kn = readline()
                        level_names[level] = kn
                    else
                        kn = level_names[level]
                    end
                else
                    if level-1 ∈ keys(level_names) && 
                       level_names[level-1] == "__key_level"
                        key_name = nothing
                    else
                        level_names[level] = "__key_level"
                    end
                    kn = key_name
                end
                other_labels[kn] = key

            else
                @warn "unhandled key_name=$key"
                other_labels["unnamed"] = key
            end
            if fields[key] !== nothing
                try
                    df = to_dataframe(fields[key]; key_name=key_name,
                                      other_labels=other_labels, kws...)
                    append!(D, df , cols=:union)
                catch
                    @warn "Hiccup"
                    @infiltrate
                end
            end
        end
        return D
    end
    function to_dataframe(fields::T where T <: NamedTuple; kws...)::DataFrame
        D = to_dataframe(Utils.namedtuple_to_dict(fields); kws...) 
        return D
    end
    function to_dataframe(X::DataFrame; other_labels=nothing, kws...)
        if other_labels !== nothing
            for (k,v) in other_labels
                X[!,k] .= v
            end
        end
        X
    end

    """
    field_to_dataframe(field::AbstractArray; other_labels=Dict())

    +purpose: converts a single field matrix into a field dataframe
    """
    function to_dataframe(F::Union{AbstractArray,Real};
            props::Vector{String}=Vector{String}([]), grid::Tuple=(),
            gridtype="center", other_labels=Dict(), name::String="value",
            explode::Bool=true, kws...)::DataFrame
        if explode
            D = ndgrid((1:size(F,i) for i in 1:ndims(F))...)
        else
            D = [1:size(F,i) for i in 1:ndims(F)]
        end
        D = OrderedDict{String,Any}("dim_$d"=>vec(D[d])
                              for d in 1:ndims(F))
        if typeof(F) <: Real
            F = [F]
        end
        if ~isempty(props)
            if grid[1] isa Tuple
                grid = Tuple([g...] for g in grid)
            end
            if gridtype == "edge"
                grid = edge_to_center.(grid)
            end
            if explode
                @info "ndgrid of grid"
                grid = ndgrid(grid...)
            end
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
