"""
    trial

trial-based transformations


outputs https://github.com/JuliaArrays/AxisArrays.jl types
"""
module tensor

    using DataFrames
    using DataStructures: OrderedDict
    using StatsBase
    using Statistics

    import Utils

    export tensor_continuous
    export quantilize, relativize

    SymStr = Union{Symbol, String}

    """
        tensor_pointproc

    create a tensor of a point process arrayed with some set of N-dimensions
    """
    function tensor_pointproc(X::DataFrame)
    end

    """
        tensor_pointproc

    create a tensor of a continuous variable from dataframe form

    ### Input
    
    `X`         -- DataFrame we'll pull a tensor out of
    `dims`      -- Columns that will become the dimension of the tensor
    `var`       -- Column that will become the measurement

    ### Output
    `T` -- Tensor; either a square array or stacked arrays of unequal dim
     
    """
    function tensor_continuous(X::DataFrame, dims::Vector, var::T where T<:SymStr;
            quantilebin::Dict{<:SymStr, <:Int}=Dict{Symbol,Int}(),
            relative::Vector{<:SymStr},
            TensType::Type=Float64)

        X = copy(X)
        
        # Clean
        X = dropmissing(X, [dims...])
        notisnan = Utils.notisnan(X[!,dims])
        notisnan = all(notisnan, dims=2)
        X = X[notisnan, :]

        X = quantilize(X, quantilebin)
        X = relativize(X, dims, relative)
        G = groupby(X[!, [dims..., var]], dims)
        counts = combine(G, nrow)
        if !(all(counts.nrow .== 1))
            TensType = Array{TensType}
            singular = false
        else
            singular = true
        end

        axes = OrderedDict(dim=>sort(unique(X[!,dim])) for dim in dims)
        sz   = [length(ax) for ax in axes]
        T = Array{TensType}(undef,sz...)
        for g in G
            searches = map((dim, ax) -> first(g)[dim] .== ax, zip(dims, axes))
            overall = (accumulate(.&, searches))[end]
            for i = 1:length(searches)
                searches[i] .&= overall
            end
            indices = Tuple(findfirst(search) for search in searches)
            T[indices...] = singular ? g[!, var] : g[!, var][1]
        end
    end

    function equalize
    end

    function quantilize(X::DataFrame, quantilebin::Dict{<:SymStr, <:Int},)
        # Rank bin anything requested
        for (var, bins) in quantilebin
            quant = X[!,var]./maximum(X[!,var])
            binned = ceil.((bins-1)*quant) 
            X[!,var] = binned
        end
    end

    function relativize(X::DataFrame, dims::Vector{<:SymStr}, var::Vector{<:SymStr}; keeporig::Bool=false)
        # Rank bin anything requested
        makerelative(x) = x .- minimum(x)
        pos = keeporig ? [v => identity => Symbol("orig"*String(v)) for v in var] : []
        combine(groupby(X, dims, sort=false), var .=> makerelative .=> var, pos...)
    end
    function relativize(X::DataFrame, dims::Vector{<:SymStr}, var::T where  T<:SymStr; keeporig::Bool=false)
        relativize(X, dims, [var]; keeporig)
    end

end
