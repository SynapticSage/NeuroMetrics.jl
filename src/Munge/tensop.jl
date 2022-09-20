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
    #axisnames

    import Utils


    function tdtw(templates::DimArray, X::DimArray; kws...)
    end
    function tdtw(templates::AbstractArray, X::DimArray, dims; kws...)
    end
    function tdtw(templates::AbstractArray, X::AbstractArray, dims; kws...)
    end

    """
        median

    get median of all values
    """
    function tmedian(X::DimArray; dims)
    end
    function tmedian(X::AbstractArray; dims)
        X = tenmat(X, row=setdiff(1:ndims(X)), col=dims)
        [Statistics.median(hcat(x...); dims=1) for x in eachrow(X)]
    end

    """
        mean

    get median of all values
    """
    function tmean(X::DimArray; dims)
    end
    function tmean(X::AbstractArray)
        X = X[:]
        X = hcat(X...)
        Statistics.mean(X; dims=2)
    end

    """
        dba

    get barycenter of all values
    """
    function dba(X::T where T<:Tensor)
    end

end
