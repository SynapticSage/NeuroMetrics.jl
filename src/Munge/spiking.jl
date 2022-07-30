module spiking

    import Munge
    using StatsBase
    import Utils: binning
    using DataFrames
    using ImageFiltering
    using AxisArrays
    using LazySets

    import Field: ReceptiveField

    export torate

    """
        nonlocality

    whether a data frame of behavior or spikes is NOT in a place field
    """
    function nonlocality(X::DataFrame, R::ReceptiveField; hull=1)::BitVector
        @assert :hullseg_grid ∈ R.metrics
        props = R.grid.props
        X = X[:, props]
        X = [element(Singleton(x)) for x in eachrow(X)]
        X .∉ VPolygon(R.metrics[:hullseg_grid][hull])
    end

    """
        locality

    whether a data frame of behavior or spikes is in a place field
    """
    function locality(X::DataFrame, R::ReceptiveField; hull=1)::BitVector
        @assert :hullseg_grid ∈ R.metrics
        props = R.grid.props
        X = X[:, props]
        X = [element(Singleton(x)) for x in eachrow(X)]
        X .∈ VPolygon(R.metrics[:hullseg_grid][hull])
    end


    function torate(spikes::DataFrame, dims=:unit; grid=nothing, kws...)
        if grid === nothing
            grid = minimum(spikes.time):0.01:maximum(spikes.time)
        elseif typeof(grid) <: Real
            grid = minimum(spikes.time):grid:maximum(spikes.time)
        end
        T =  Munge.tensorize(spikes, [dims], :time)
        M = [torate(x; grid, kws...) for x in T]
        timeax = M[1].axes
        M = hcat(M...)
        neuronax = T.axes
        AxisArray(M, timeax..., neuronax...)
    end

    function torate_windowdia(times::AbstractArray; grid, windowsize::Real,
            gaussian::Real=0)
        torate_windowrad(times; grid, radius=windowsize/2, gaussian)
    end
    #function torate(times::AbstractArray; grid, windowsizes::Tuple{<:Real,<:Real})
    #end
    function torate_windowrad(times::AbstractArray; grid, radius::Real, 
                           gaussian::Real=0)
        count = binning.inside(times, grid, radius)
        ker = Kernel.gaussian((gaussian,))
        δ  = grid[2] - grid[1]
        centers = collect(grid)
        AxisArray(imfilter(count ./ δ, ker),
                  Axis{:time}(centers))
    end
    function torate(times::AbstractArray; grid, gaussian=0.25)::AxisArray
        δ  = grid[2] - grid[1]
        count = fit(Histogram, vec(times), grid)
        gaussian = gaussian * 0.1/δ
        ker = Kernel.gaussian((gaussian,))
        centers = binning.edge_to_center(collect(grid))
        AxisArray(imfilter(count.weights ./ δ, ker),
                  Axis{:time}(centers))
    end

    # CELL COFIRING
    function xcorr(rate::AxisArray, cell1::T , cell2::T; lags=-20:20) where 
        T<:Union{Int16,Int32,Int64}
        x, y = rate[Axis{:unit}(cell1)],
               rate[Axis{:unit}(cell2)]
        StatsBase.crosscor(x, y, lags)
    end

    function xcorr(spikes::DataFrame; lags=-200:200, kws...)
        units1 = unique(spikes.unit)
        units2 = unique(spikes.unit)
        R = torate(spikes; kws...)
        results = []
        for (cell1,cell2) in Iterators.product(units1,units2)
            push!(results, xcorr(R, cell1, cell2; lags))
        end
    end

end

