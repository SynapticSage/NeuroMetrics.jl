#=
Methods for dealing with keys that summarize
batches of timeshift measurements under different conditions

Keys of these Dict objects tend to contain NamedTuples that describe
the conditions of the measurement
=#

"""
For converting a single key
"""
function keytostring(key::NamedTuple)
    key = Dict(zip(keys(key), values(key)))
    for (k,v) in key
        tv = typeof(v)
        if tv  <: Vector{AbstractString}
            v = join(tv, "-")
        elseif tv <: Real
            v = round(tv; sigdigits=2)
        end
        key[k] = v
    end
    string(key)
end

"""
For converting whole dicitionaries containing these keys
"""
function keytostring(D::Dict{<:NamedTuple, <:Any})
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
