module operation
using Statistics
using NaNStatistics
using Bootstrap
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
                     isdensity::Bool=false)
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
        end
    end
    return field
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

end
