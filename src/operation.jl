
module operation
using Statistics
using NaNStatistics
using Bootstrap
include("utils.jl")
export utils

"""
field_operation

used to collapse a field onto 1D
"""
function unary(field, operation::Function=mean;
    kws...)
    if typeof(field) <: AbstractArray
        field = operation(field; kws...)
    elseif field isa Dict
        field = Dict(name=>unary(f, operation; kws...) for (name,f) in
                     field)
    end
    return field
end
"""
field_operation

used to do an operation between fields
"""
function binary(fieldA, fieldB, operation::Function=(./); kws...)
    if typeof(fieldA) <: AbstractArray
        field = operation(fieldA, fieldB; kws...)
    elseif fieldA isa Dict
        field = Dict(namesA=>binary(fieldA[namesA], fieldB[namesA],
                                  operation; kws...) for namesA in
                     keys(fieldA))
    else
        throw(TypeError("fielld is wrong type $(typeof(fieldA))"))
    end
    return field
end

isdens = Dict(:hist=>false, :kde=>false, :beh=>false, :behdens=>true)

"""
marginalize(field, dims)

example: p(x,y,Γ) ->(Γ) p(x,y)
"""
function marginalize(field::NamedTuple; dims::Vector{Int}=[],
                     isdensity::Bool=false, dosqueeze::Bool=false)
    field = Dict(pairs(field))
    for (key, item) ∈ field
        if item == nothing
            continue
        end
        if occursin("grid", String(key))
            S = setdiff(1:length(item), dims)
            item = item[S]
        else
            dens = isdensity ? true : isdens[key]
            item = marginalize(item; dims=dims, isdensity=dens)
            if dosqueeze; item = squeeze(item); end
        end
        field[key] = item
    end
    return (;field...)
end
function marginalize(fields::Dict; dims::Vector{Int}=[],
        isdensity::Bool=false)
    return Dict(key=>marginalize(field, dims=dims, isdensity=isdensity)
                for (key,field) in fields)
end

"""
marginalize

TODO renan areas where all() missing values in relevent dims
"""
function marginalize(field::AbstractArray; dims::Vector{Int}=[],
        isdensity::Bool=false)
    for dim in dims
        field = nansum(field, dims=dim)
    end
    if isdensity; field = field./nansum(field); end
    return field
end


"""
squeeze

meethods for squeezing field objects and field dicts
"""
function squeeze(field::AbstractArray)
    return utils.squeeze(field)
end
function squeeze(fields::Dict)
    return Dict(key=>utils.squeeze(field) for (key,field) in fields)
end
function squeeze(field::NamedTuple)
    field = Dict(pairs(field))
    for (key, item) ∈ field
        if item == nothing
            continue
        end
        if !(occursin("grid", String(key)))
            item = squeeze(item)
        end
        field[key] = item
    end
    return (;field...)
end

function apply(func::Function, D::Dict...; keypolicy::Function=intersect,
        kws...)
    K = [Set(keys(d)) for d∈D]
    K = accumulate(keypolicy, K)[end]
    result = Dict()
    for k∈K
        d = (d(k) for d in D)
        result[k] = func(d...; kws...)
    end
    return result
end
function apply(func::Function, D::Dict, X::AbstractArray; kws...)
    result = Dict()
    for k∈keys(D)
        result[k] = func(D, X; kws...)
    end
    return result
end
function apply(func::Function, X::AbstractArray, D::Dict; kws...)
    result = Dict()
    for k∈keys(D)
        result[k] = func(X, D[k]; kws...)
    end
    return result
end

end
