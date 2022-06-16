module operation

using Statistics
using NaNStatistics
using Bootstrap
using DataStructures
using Utils

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
        @info operation
        field = operation(fieldA, fieldB; kws...)
    elseif typeof(fieldA) <: Real
        @info operation
        field = operation(fieldA, fieldB; kws...)
        @info field
    elseif typeof(fieldA) <: AbstractDict
        K = intersect(keys(fieldA), keys(fieldB))
        field = Dict(namesA=>binary(fieldA[namesA], fieldB[namesA],
                                  operation; kws...) for namesA in K)
    else
        throw(TypeError("operation::binary", Union{AbstractArray, AbstractDict}, typeof(fieldA)))
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
function marginalize(F::NamedTuple; dims::Vector{Int}=[],
                     isdensity::Bool=false, dosqueeze::Bool=true,
                     removecount::Bool=true, debug::Bool=false)
    F = Dict(pairs(F))
    for (key, item) ∈ F
        if item == nothing
            continue
        end
        #println(key)
        sKey = String(key)
        if endswith(sKey, "sq")
            continue
        end
        println(sKey)
        if item == nothing
           continue 
        elseif occursin("grid", sKey)
            S = setdiff(1:length(item), dims)
            item = item[S]
        elseif startswith(sKey, "R") # rates
            dens = isdensity ? true : isdens[key]
            item = marginalize(item; dims=dims, isdensity=dens,  #TODO should item be count Cₓ instead of rate Rₓ
                            normalizeFR=F[:occ])
        elseif (startswith(sKey, "C") && !(removecount)) || sKey == "occ" 
            dens = isdensity ? true : isdens[key]
            item = marginalize(item; dims=dims, isdensity=dens)
        elseif sKey == "dims"
            item = item[setdiff(1:length(item), dims)]
        end
        F[key] = item
    end
    if removecount
        if debug; println("Removing counts"); end
        for key ∈ keys(F)
            if startswith(String(key), "C")
                delete!(F, key)
            end
        end
    end
    if dosqueeze 
        if debug; println("Squeezing"); end
        for key ∈ String.(keys(F))
            if endswith(key, "sq")
                continue
            end
            if any(startswith.(key,["R","C"]))
                item = F[Symbol(key)]
                if debug; println("Size item=$(size(si(item)))"); end
                item = operation.squeeze(item);
                if debug; println("Size item=$(size(si(item)))"); end
                F[Symbol(key * "sq")] = item
            end
        end
    end
    return NamedTuple(F)
end
function marginalize(fields::AbstractDict; kws...)
    T = typeof(fields)
    return T(key=>marginalize(field; kws...) for (key,field) in fields)
end
function marginalize(field::AbstractArray; dims::Vector{Int}=[],
        isdensity::Bool=false, 
        normalizeFR::Union{Nothing, AbstractArray}=nothing,
        maintainNaN::Bool=true, debug::Bool=false)
    donorm = !(normalizeFR isa Nothing)
    if debug; println("(1) size(field) = $(size(field)), dims=$dims"); end
    if donorm
        occzero = normalizeFR .== NaN
        field[occzero] .= NaN
    end
    for dim in dims
        field = nansum(field, dims=dim)
        if donorm
            normalizeFR = nansum(normalizeFR, dims=dim)
        end
        if debug; println("(2) dim=$dim => size(field) = $(size(field))"); end
    end
    if donorm
        field = field ./ normalizeFR
    end
    if isdensity; field = field./nansum(field); end
    if debug; println("(3) size(field) = $(size(field))"); end
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
            item = operation.squeeze(item)
        end
        field[key] = item
    end
    return (;field...)
end

"""
Apply an operation to a dict or set of dict objects containing some type of
field information
"""
function apply(func::Function, D::Dict...; tolerate_error::Bool=false,
        keypolicy::Function=intersect, kws...)
    #println("apply 3")
    K = [Set(keys(d)) for d∈D]
    K = accumulate(keypolicy, K)[end]
    result = Dict()
    for k∈K
        d = (d[k] for d in D)
        if tolerate_error
            try
                result[k] = func(d...; kws...)
            catch
                @warn "key=$k during operation.apply results in error"
                result[k] = nothing
            end
        else
            result[k] = func(d...; kws...)
        end
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
    if data[1] isa Int64
        X = convert(Array{Float64}, data);
    elseif data[1] isa Int32
        X = convert(Array{Float32}, data);
    else
        X = data
    end
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

function cast32(F::NamedTuple; kws...)
    F = OrderedDict(pairs(F))
    for d in keys(F)
        if F[d] == nothing
            continue
        end
        println(d)
        if F[d] isa Dict
            F[d] = operation.apply(x->Float32.(x), F[d]; kws...)
        elseif F[d] isa AbstractArray{Float64}
            F[d] = Float32.(F[d])
        elseif F[d] isa AbstractArray{Int64}
            F[d] = Int32.(F[d])
        end
    end
    return NamedTuple(F)
end
function selectrand(D::AbstractDict, N::Int)
    K = Tuple(keys(D))
    T = typeof(D)
    r = rand(1:length(K), N)
    K = K[r]
    D  = T(k=>D[k] for k in K)
end
const sr = selectrand

function selectkey(D::AbstractDict, P::Pair...)
    T = typeof(D)
    K = keys(D)
    for (tupfield, command) in P
        if tupfield isa String
            tupfield = Symbol(tupfield)
        end
        if command isa Function
           if tupfield isa Sybmol
                func = command
                K = filter(k->func(k[tupfield]), K)
           elseif tupfield == All()
               K = filter(func, K)
           end
        else
            obj = command
            K = filter(k->k[tupfield] == obj, K)
        end
    end
    T(k=>D[k] for k in K)
end
const sk = selectkey

function selectind(D::Dict, ind::Int=1)
    K = Tuple(keys(D))
    D[K[ind]]
end
const si = selectind

function split_dimofval(F::Dict, dim::Int; 
        name::Union{String, Symbol}=:auto,
        indices::Union{Vector,Tuple,Nothing}=nothing)
    F′ = Dict()
    if name isa Vector{String}
        name = name[dim] # user passed in fields.dims
        if indices != nothing
            indices = indices[dim] # user passesd in fields.grid
        end
    end
    if name isa String; name=Symbol(name); end
    for (key, value) in F
        @assert key isa NamedTuple
        f = F[key]
        if indices == nothing
            indices = 1:size(value,dim)
        end
        slices = eachslice(value, dims=dim)
        for (index, slice) in zip(indices, slices)
            key′ = (; key..., Dict(name=>index)...)
            F′[key′] = slice
        end
    end
    F′
end

"""
group

all field objects who share a keey are grouped by their parent fields

the keys for the parent fields are passedd in `names`
"""
function group(fields::Union{Tuple,Vector},
        names::Union{Vector{String},Vector{Symbol}};
        grouptype::Type=Dict, requireAll=false)

    keyset = keys(fields[1])
    for i in range(2,length(fields))
        keyset = union(keyset, keys(fields[i]))
    end
    keyset = Tuple(keyset)
    Grouping = Dict{typeof(keyset[1]), Any}()
    for key in keyset
        if requireAll
            notin = [!(key in keys(field)) for (name,field) in
                     zip(names,fields)]
            if any(notin)
                continue
            end
        end
        if grouptype <: AbstractArray
            Grouping[key] = Dict(name=>field[key]
                                  for (name,field) in zip(names, fields)
                                  if key in keys(field))
        elseif grouptype == NamedTuple
            names = Symbol.(names)
            Grouping[key] = (;((name,field[key]) for (name,field) in
                               zip(names, fields) if key in keys(field))...)
        else
            Grouping[key] = (grouptype)(field[key]
                                  for (name,field) in zip(names, fields)
                                  if key in keys(field))
        end
    end
    return Grouping
end

function catfields
end

function reorderfields
end


end
