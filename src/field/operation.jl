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

isdens = Dict(:Cₕ=>false, :Cₖ=>false, :occ=>false, :behdens=>true,
              :occR=>true, :Rₕ=>false, :Rₖ=>false)

"""
marginalize(field, dims)
example: p(x,y,Γ) ->(Γ) p(x,y)

marginalizes a field object (result of get_fields)

or a dict of fields

or a field itself (AbstractArray)
"""
function marginalize(field::NamedTuple; dims::Vector{Int}=[],
                     isdensity::Bool=false, dosqueeze::Bool=false)
    field = Dict(pairs(field))
    for (key, item) ∈ field
        if item == nothing
            continue
        end
        println(key)
        sKey = String(key)
        if occursin("grid", sKey)
            S = setdiff(1:length(item), dims)
            item = item[S]
        elseif startswith(sKey, "R") # rates
            dens = isdensity ? true : isdens[key]
            item = marginalize(item; dims=dims, isdensity=dens, 
                            normalizeFR=field[:occ])
            if dosqueeze; item = squeeze(item); end
        elseif startswith(sKey, "C") || sKey == "occ"
            dens = isdensity ? true : isdens[key]
            item = marginalize(item; dims=dims, isdensity=dens)
            if dosqueeze && sKey != "occ" ; item = squeeze(item); end
        end
        if dosqueeze && "occ" in keys(field)
        end
        field[key] = item
    end
    return NamedTuple(field)
end
function marginalize(fields::Dict; kws...)
    return Dict(key=>marginalize(field; kws...) for (key,field) in fields)
end
"""
marginalize

TODO renan areas where all() missing values in relevent dims
"""
function marginalize(field::AbstractArray; dims::Vector{Int}=[],
        isdensity::Bool=false, 
        normalizeFR::Union{Nothing, AbstractArray}=nothing,
        maintainNaN::Bool=true)
    for dim in dims
        field = nansum(field, dims=dim)
    end
    if !(normalizeFR isa Nothing)
        occzero = normalizeFR .== NaN
        field = field ./ normalizeFR
        field[occzero] .= NaN
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

"""
Apply an operation to a dict or set of dict objects containing some type of
field information
"""
function apply(func::Function, D::Dict...; keypolicy::Function=intersect,
        kws...)
    #println("apply 3")
    K = [Set(keys(d)) for d∈D]
    K = accumulate(keypolicy, K)[end]
    result = Dict()
    for k∈K
        d = (d[k] for d in D)
        result[k] = func(d...; kws...)
    end
    return result
end
function apply(func::Function, D::Dict, X::AbstractArray; kws...)
    #println("apply 1")
    result = Dict()
    for k∈keys(D)
        result[k] = func(D[k], X; kws...)
    end
    return result
end
function apply(func::Function, X::AbstractArray, D::Dict; kws...)
    #println("apply 2")
    result = Dict()
    for k∈keys(D)
        result[k] = func(X, D[k]; kws...)
    end
    return result
end

function occnorm(F::NamedTuple, dfields=["Cₕ","Cₖ"]; occfield="occ",
        ozifield="occzeroinds", kws...)
    F = Dict(pairs(F))
    for d in dfields
        newName = replace(d, "C"=>"R")
        isANormRateField = startswith(d, "R")
        noData = F[Symbol(d)] isa Nothing 
        if noData || isANormRateField
            continue
        end
        F[Symbol(newName)] = operation.apply(operation.occnorm, 
                                             F[Symbol(d)],
                                             F[Symbol(occfield)]; 
                                             occzeroinds=F[Symbol(ozifield)],
                                             kws...)
    end
    return NamedTuple(F)
end
function occnorm(data::AbstractArray,
        occupancy::Union{Array,Matrix};
        occzeroinds::Union{AbstractArray,Nothing}=nothing,
        gaussian::Real=0.0)
    X = convert(Array{Float64}, data);
    X = replace(X, NaN=>0)
    X = X ./ occupancy;
    if gaussian != 0.0
        X = imfilter(X, Kernel.gaussian(gaussian))
    end
    if occzeroinds != nothing
        #println("apply occzeroinds")
        X[occzeroinds] .= NaN;
    end
    return X
end

end
