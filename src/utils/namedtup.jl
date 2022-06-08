
#=
This sub library is about functions for named
tuples and for dicts who have named tuple
keys.
=#

function namedtuple_to_dict(X::NamedTuple)
    Dict(zip(keys(X), values(X)))
end
function pop(X::NamedTuple, key)
    X = namedtuple_to_dict(X)
    pop!(X, key)
    NamedTuple(X)
end
function namedtupkeys_to_df(K::Base.KeySet)
    df = DataFrame()
    for k in K
        k = Dict(zip(keys(k), values(k)))
        k = DataFrame(k)
        append!(df, k)
    end
    df
end

function remove_key_item(k::NamedTuple, item)
    k = Dict(zip(keys(k), values(k)))
    k = remove_key_item(k, item)
    NamedTuple(k)
end

# Dicts with named tuple keys
# May eventually split this off to its own sublib

function namedtupkeys_to_df(D::AbstractDict)
    K = keys(D)
    namedtupkeys_to_df(K)
end

"""
Counts the number of matches a dict with named tuple keys to a search namedtuple
"""
function countmatch(D::AbstractDict{<:NamedTuple,Any}, 
                             search::NamedTuple)
    countmatch(keys(D), search)
end

"""
Counts the number of matches in vector of namedTuples to a search namedtuple
"""
function countmatch(N::Vector{<:NamedTuple}, search::NamedTuple)
    results = []
    for n in N
        key_results = []
        keys_intersect = intersect(keys(n), keys(search))
        for k in keys_intersect
            match = n[k] == search[k]
            push!(key_results, match)
        end
        push!(results,sum(key_results))
    end
    return results
end
function bestpartialmatch(K::Vector{<:NamedTuple}, search::NamedTuple)
    argmax(countmatch(K, search))
end
function bestpartialmatch(D::AbstractDict{<:NamedTuple, Any}, search:NamedTuple)
    argmax(countmatch(D, search))
end
