module Keys

    using ProgressMeter
    using DataFrames
    using Infiltrator
    using DataStructures: OrderedDict
    import DIutils
    export pluto_executekey, pluto_keyselector
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

    function pluto_keyselector(data::Dict)
        # Determine which properties CAN be toggled
        allsets   = collect(keys(data))
        totalkeys = union(keys.(allsets)...)
        uvals = OrderedDict()
        for key in totalkeys
            U = unique([getindex(s, key) for s in allsets
                        if key ∈ propertynames(s)])
            push!(uvals, key=>U)
        end

        # Determine how they would be represented as controls and actual objects
        actuals  = OrderedDict()
        controls = OrderedDict()
        rep(x) = replace(x, "\""=>"", "]"=>"", "["=>"")
        for (K,V) in uvals
            sortV = [sort(V)..., nothing]
            push!(controls, K => [rep("$v") for v in sortV])
            for v in sortV
                push!(actuals,  ((v === nothing) ? "nothing" : "$v") => v)
            end
        end
        controls, actuals
    end

    function pluto_executekey(data::AbstractDict, actuals::AbstractDict, specification::Pair...;
            hard_settings=(;grid=:adaptive, first=-2.0, step=0.05, last=2.0), return_findthis=false)

        soft_settings = NamedTuple(Dict(k=>actuals[v] for (k,v) in specification))
        
        find_this = (;hard_settings..., soft_settings...)
        key=DIutils.namedtup.bestpartialmatch(keys(data), find_this;
                                  nothing_means_removekey=true)
        return_findthis ? (key, find_this) : key
    end

end
