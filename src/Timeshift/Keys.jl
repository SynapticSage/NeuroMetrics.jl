module Keys

    using ProgressMeter
    using DataFrames
    using Infiltrator
    #=
    Methods for dealing with keys that summarize
    batches of timeshift measurements under different conditions

    Keys of these Dict objects tend to contain NamedTuples that describe
    the conditions of the measurement
    =#

    #=
    ------------------
    ------------------
    KEYS TO STRING KEYS
    ------------------
    ------------------
    =#

    """
    For converting a single key
    """
    function keytostring(key::NamedTuple, extra_replacements::Pair...; remove_numerics::Bool=false)
        key = Dict(zip(keys(key), values(key)))
        for (k,v) in key
            tv = typeof(v)
            if tv  <: Vector{<:AbstractString}
                v = join(v, "-")
                key[k] = v
            elseif tv <: Real
                if remove_numerics
                    pop!(key, k)
                else
                    v = round(v; sigdigits=2)
                    key[k] = v
                end
            else
                key[k] = v
            end
        end
        replace(string(NamedTuple(key)), " = "=>"=", "\""=>"", ")"=>"", "("=>"",
               extra_replacements...)
    end

    """
    For converting whole dicitionaries containing these keys
    """
    function keytostring(D::Dict{<:NamedTuple, <:Any}, extra_replacements::Pair...; kws...)
        Dict(keytostring(k, extra_replacements...; kws...)=>v for (k,v) in D)
    end

    """
    Even shorter than the keytostring method
    """
    function shortcutstring(key::NamedTuple; remove_numerics::Bool=true)
        val_list = []
        key = Dict(zip(keys(key),values(key)))
        for (k,v) in key
            tv = typeof(v)
            if tv  <: Vector{<:AbstractString}
                v = join(v, "-")
            elseif tv <: Real
                if !(remove_numerics)
                    v = round(v; sigdigits=2)
                    v = string(v)
                else
                    v = ""
                end
            elseif tv === Symbol
                v = ":" * string(v)
            else
                v=string(v)
            end
            push!(val_list, string(v))
        end
        replace(join(val_list,""), " = "=>"=", "\""=>"")
    end

    #=
    ONE TIME CORRECTIONS

    Corrections for key structures that were at one point different
    =#

    function _correct_mains(D::Dict{<:NamedTuple, <:Any})
        kD = keys(D)
        @showprogress for (i,nt_key) in enumerate(kD)
            @info "key=$key \n i=$i"
            @infiltrate
            value = pop!(D, nt_key)
            nt_key = Dict(zip(keys(nt_key), values(nt_key)))
            if :shuf ∈ keys(nt_key)
                pop!(nt_key, :shuf)
            end
            if :shuffle_type ∈ keys(nt_key)
                pop!(nt_key, :shuffle_type)
            end
            @infiltrate
            new_key = NamedTuple(nt_key)
            D[new_key] = value
        end
    end

end
