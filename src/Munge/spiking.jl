module spiking

    import Munge
    using StatsBase
    import Utils: binning
    using DataFrames
    using ImageFiltering
    using AxisArrays

    import Field: ReceptiveField

    function nonlocality(times::Vector, beh::DataFrame, R::ReceptiveField)
    end
    function nonlocality(beh_of_spikes::Array, R::ReceptiveField)
    end

    function torate(spikes::DataFrame, dims=:unit; kws...)
        T =  Munge.tensorize(spikes, [dims], :time)
        M = [torate(x; kws...) for x in T]
        timeax = M[1].axes
        M = hcat(M...)
        neuronax = T.axes
        AxisArray(M, timeax..., neuronax...)
    end
    function torate(times::AbstractArray; grid, windowsize::Tuple{<:Real,<:Real})
    end
    function torate(times::AbstractArray; grid, windowsize::Real, gaussian::Real)
        torate(times; grid, radius=windowsize/2)
    end
    function torate(times::AbstractArray; grid, radius::Real)
        binning.inside(times, grid, radius)
    end
    function torate(times::AbstractArray; grid, gaussian=0.25)::AxisArray
        δ  = grid[2]-grid[1]
        count = fit(Histogram, vec(times), grid)
        gaussian = gaussian * 0.1/δ
        ker = Kernel.gaussian((gaussian,))
        centers = binning.edge_to_center(collect(grid))
        AxisArray(imfilter(count.weights ./ δ, ker),
                  Axis{:time}(centers))
    end



end

