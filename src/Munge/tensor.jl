"""
    trial

trial-based transformations


outputs https://github.com/JuliaArrays/DimensionalData.jl types
"""
module tensor

    using DataFrames
    using DataStructures: OrderedDict
    using StatsBase
    using Statistics
    import Statistics: mean, median
    using Infiltrator
    using ProgressMeter
    using Missings
    using DimensionalData
    using DimensionalData.Dimensions
    using Interpolations
    using DynamicAxisWarping
    using TensorToolbox
    using CategoricalArrays

    import DIutils

    export tensorize, tensorize_continuous
    export quantilize, relativize, equalize, gravity, digitize
    export tenmat

    SymStr          = Union{Symbol, String, Int} # allowed column names
    Tensor          = DimArray{<:Union{Missing,Array,Real}}
    DataFrameTensor = DimArray{<:Union{Missing,DataFrame}}
    AllTensor       = Union{Tensor, DataFrameTensor}


    """
        tensorize

    create a tensor of a continuous variable from dataframe form

    ### Input
    
    `X`         -- DataFrame we'll pull a tensor out of
    `dims`      -- Columns that will become the dimension of the tensor
    `var`       -- Column that will become the measurement

    ### Output
    `T` -- Tensor; either a square array or stacked arrays of unequal dim
     
    """
    function tensorize(X::DataFrame, dims::Vector{T} where T<:SymStr, 
            var::Union{T, Vector{T}} where T<:SymStr;
            TensType::Type=Float64)::AllTensor

        if !(var isa Vector)
            var = [var]
        end

        X = copy(X)
        X = clean(X, dims)

        group_var = union(dims, var)
        G = groupby(X[!, group_var], dims)
        counts = combine(G, nrow)
        if !(all(counts.nrow .== 1))
            TensType = Array{TensType}
            scalar = false
        else
            scalar = true
        end
        arraycols = any(eltype.(eachcol(X[!,var])) .<: Array) ?
                    true : false

        axs = OrderedDict(dim=>sort(unique(X[!,dim])) for dim in dims)
        sz   = [length(ax) for (_,ax) in axs]
        TensType = TensType !== DataFrame ? 
           Array{Union{Missing,TensType}} :
           Union{Missing,DataFrame}

        T = TensType(missing,sz...)

        @showprogress for g in G
            # Lookup the dim values of g in the axs lists
            function match_func(dim, ax)
                first(g)[dim] .== ax
            end
            searches = map(match_func, keys(axs),values(axs))
            #overall = accumulate(.&, searches)[end]
            #for i = 1:length(searches)
            #    searches[i] .&= overall
            #end
            indices = findfirst.(searches)
            val = if scalar 
                g[!, var][1]
            elseif arraycols
                G = []
                for col in eachcol(g[!,var])
                    if eltype(col) <: Array
                        [push!(G,c) for c in eachcol(hcat(col...)')]
                    else
                        push!(G,col) 
                    end
                end
                hcat(G...)
            else
                Array(g[!,var])
            end
            T[indices...] = val
        end
        DimArray(T, Tuple(Dim{Symbol(name)}(index isa String ? categorical(index) : index)
                          for (name,index) in axs))
    end

    """
        interp

    when a tensor is a tensor{array}, equalize the number of samples across
    all by some kind of interpolation (using Interpolate arguments from
    Interpolations.jl)
    """
    function interp(X::AbstractArray{<:Array}, interp=nothing, 
                extrap=OnGrid())::AbstractArray
        if interp === nothing
            interp =  ndims(X[1]) == 1 ? 
                BSpline(Linear()) : (BSpline(Linear()), NoInterp())
        end
        szs = [size(t,1) for t in X]
        max_sz = maximum(szs)
        for (i, x) in zip(eachindex(X), X)
            @infiltrate
            int = interpolate(x, interp, extrap)
            t = int[range(1, size(x,1), length=max_sz),:]
            X[i] = t
        end
        X
    end
    

    function quantilize(X::DataFrame, dims::Vector{<:SymStr},
            quantilebin::Dict{<:SymStr, <:Int},)::DataFrame
        X = clean(X, dims)
        # Rank bin anything requested
        for (var, bins) in quantilebin
            quant = X[!,var]./maximum(X[!,var])
            binned = ceil.((bins-1)*quant) 
            X[!,var] = binned
        end
        X
    end

    """
    Digitize columns
    """
    function digitize(X::DataFrame, dims::Vector{<:SymStr}, bins::Union{Vector{<:Int},Int})::DataFrame
        bins = bins isa Int ? fill(bins, size(dims)) : bins
        X = clean(X, dims)
        # Rank bin anything requested
        for (dim,bin) in zip(dims, bins)
            X[!,dim] = DIutils.binning.digitize(X[!,dim], bin)
        end
        X
    end

    """
    relativize
    
    makes a measurement relative to a center or  minimum within 
    each group (dims)
    """
    function relativize(X::DataFrame, groupvars::Vector{<:SymStr},
            var::Vector{<:SymStr}; keeporig::Bool=false)::DataFrame
        X = clean(X, groupvars)
        # Rank bin anything requested
        makerelative(x) = x .- minimum(x)
        pos = keeporig ? [v => identity => Symbol("orig"*String(v)) for v in var] : []
        combine(groupby(X, groupvars, sort=false), var .=> makerelative .=> var, pos...)
    end
    function relativize(X::DataFrame, groupvars::Vector{<:SymStr}, var::T where
            T<:SymStr; keeporig::Bool=false)::DataFrame
        relativize(X, groupvars, [var]; keeporig)
    end

    function clean(X::DataFrame, dims::Vector{<:SymStr}=All())::DataFrame
        # Clean
        X = dropmissing(X, [dims...])
        notisnan = hcat(DIutils.notisnan.(eachcol(X[!,dims]))...)
        notisnan = DIutils.squeeze(all(notisnan, dims=2))
        X = X[notisnan, :]
    end


    #                                          
    #            --.--                         
    #              |  ,---.,---.,---.,---.,---.
    #              |  |---'|   |`---.|   ||    
    #              `  `---'`   '`---'`---'`    
    #                                          
    #                                                     
    #                           |              o          
    #            ,---.,---.,---.|---.,---.,---..,---.,---.
    #            |    |---'`---.|   |,---||   |||   ||   |
    #            `    `---'`---'`   '`---^|---'``   '`---|
    #                                     |          `---'

    # --- Overloading TensorToolbox.jl functions to work with array of array
    # BORROWED FROM TensorToolbox.jl
        """
        tenmat(X,n)
        tenmat(X,row=[],col=[])
        tenmat(X::ttensor,n)
        tenmat(X::ktensor,n)
    Mode-n matricization of a tensor or matricization by row and column vectors R and C.
    """
    function TensorToolbox.tenmat(X::AbstractArray{<:Union{AbstractArray,Missing},N}, n::Integer) where N
        @assert(n<=ndims(X),"Mode exceedes number of dimensions")
        sz=size(X)
        m=setdiff(1:N,n)
        if [n;m]!=collect(1:N)
            X=permutedims(X,[n;m])
        end
        reshape(X,sz[n],prod(sz[m]))
    end

    function TensorToolbox.tenmat(X::AbstractArray{<:Union{AbstractArray,Missing},N};row=[],col=[]) where N
        @assert(row!=[] || col!=[],"Al least one of row and col needs to be specified.")
        if row!=[] && col!=[]
            @assert(sort([row;col])==collect(1:N),"Incorrect mode partitioning.")
        elseif row==[]
            @assert(!(false in [c in collect(1:N) for c in col]),"Incorrect modes.")
            if isa(col,Integer)
                col=[col]
            end
            row=collect(1:N)
            deleteat!(row,sort(col))
        else
            @assert(!(false in [r in collect(1:N) for r in row]),"Incorrect modes.")
            if isa(row,Integer)
                row=[row]
            end
            col=collect(1:N)
            deleteat!(col,sort(row))
        end
        sz=size(X)
        J=prod(sz[row]);K=prod(sz[col])
        if [row;col]!=collect(1:N)
            X=permutedims(X,[row;col])
        end
        reshape(X,J,K)
    end

    """
    matten(A,n,dims)
    matten(A,R,C,dims)
    Fold matrix A into a tensor of dimension dims by mode n or by row and column vectors R and C.
    """
    function matten(A::AbstractMatrix{<:Union{AbstractArray,Missing}},
            n::Integer, dims::AbstractVector{<:Integer})
        @assert(dims[n]==size(A,1),"Dimensions mismatch")
        m = setdiff(1:length(dims), n)
        @assert prod(dims[m])==size(A,2)
        X = reshape(A,[dims[n];dims[m]]...)
        permutedims(X,invperm([n;m]))
    end
    function matten(A::AbstractMatrix{<:Union{AbstractArray,Missing}},
            row::Vector{<:Integer}, col::AbstractVector{<:Integer},
            dims::Vector{<:Integer})
        @assert(prod(dims[row])==size(A,1) && prod(dims[col])==size(A,2),"Dimensions mismatch")
        X = reshape(A,[dims[row];dims[col]]...)
        permutedims(X,invperm([row;col]))
    end

    #"""
    #    equalize

    #equalizes dimensions -- the block will have missing values where no samples
    #are had. this attempts to make a core dense block with no missing values by
    #throwing out values from slices, so all slices have equal samples and no
    #missing values.
    #"""
    #function equalize(x::AbstractArray, dim; missval=missing, thresh=minimum)

    #    """
    #    Cause all of the non-missing values to "fall" to
    #    the beginning of the sequence in their proper order
    #    """
    #    function gravity(y)
    #        new_y = similar(y)
    #        miss = y .=== missval
    #        good = (!).(miss)
    #        new_y[1:sum(good)] = y[good]
    #        new_y[sum(good)+1:length(y)] = y[miss]
    #        y
    #    end

    #    szx = [size(x)...]
    #    x = tenmat(x, dim)
    #    target_num = thresh([sum(s .=== missval) for s in eachslice(x; dims=2)])
    #    x  = hcat(gravity.(eachcol(x))...)
    #    x = x[1:target_num, :]
    #    szx[dim] = target_num
    #    return matten(x, dim, szx)
    #    
    #end

    #"""
    #    equalize

    #equalizes to create an ND block of nonmissing values. this creates
    #a balanced number of samples across all axes
    #"""
    #function equalize(x::DimArray, dim; missval=missing, thresh=minimum)
    #    axs = x.dims
    #    if !(dim isa Int)
    #        dim = findfirst(dim .== name.(DimensionalData.dims(x)))
    #    end
    #    data = equalize(x.data, dim; missval, thresh)
    #    axs = Tuple(i != dim ? axs[i] : 
    #                Dim{name(DimensionalData.dims(axs[i]))}(axs[i][1:size(data,dim)]) 
    #           for i in 1:length(axs))
    #    isempty(data) ? @error("empty data") : nothing
    #    DimArray(data, axs)
    #end


end
