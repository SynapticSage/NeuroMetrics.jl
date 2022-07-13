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
    using Infiltrator
    using ProgressMeter
    using Missings

    import Utils

    export tensor_continuous
    export quantilize, relativize

    SymStr = Union{Symbol, String}

    """
        tensor_pointproc

    create a tensor of a point process arrayed with some set of N-dimensions
    """
    function tensor_pointproc(X::DataFrame)::Array
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
            TensType::Type=Float64)::Array

        X = copy(X)
        X = clean(X, dims)

        G = groupby(X[!, [dims..., var]], dims)
        counts = combine(G, nrow)
        if !(all(counts.nrow .== 1))
            TensType = Array{TensType}
            singular = false
        else
            singular = true
        end

        axs = OrderedDict(dim=>sort(unique(X[!,dim])) for dim in dims)
        sz   = [length(ax) for (_,ax) in axs]
        TensType = TensType !== DataFrame ? 
           Array{Union{Missing,TensType}}(missing,sz...) :
           Union{Missing,DataFrame}

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
            T[indices...] = singular ? g[!, var][1] : Array(g[!, var])
        end
        T
    end

    """
        equalize

    equalizes dimensions -- the block will have missing values where no samples
    are had. this attempts to make a core dense block with no missing values by
    throwing out values from slices, so all slices have equal samples and no
    missing values.
    """
    function equalize
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

    function relativize(X::DataFrame, dims::Vector{<:SymStr},
            var::Vector{<:SymStr}; keeporig::Bool=false)::DataFrame
        X = clean(X, dims)
        # Rank bin anything requested
        makerelative(x) = x .- minimum(x)
        pos = keeporig ? [v => identity => Symbol("orig"*String(v)) for v in var] : []
        combine(groupby(X, dims, sort=false), var .=> makerelative .=> var, pos...)
    end
    function relativize(X::DataFrame, dims::Vector{<:SymStr}, var::T where
            T<:SymStr; keeporig::Bool=false)::DataFrame
        relativize(X, dims, [var]; keeporig)
    end

    function clean(X::DataFrame, dims::Vector{<:SymStr}=All())::DataFrame
        # Clean
        X = dropmissing(X, [dims...])
        notisnan = hcat(Utils.notisnan.(eachcol(X[!,dims]))...)
        notisnan = Utils.squeeze(all(notisnan, dims=2))
        X = X[notisnan, :]
    end

end
