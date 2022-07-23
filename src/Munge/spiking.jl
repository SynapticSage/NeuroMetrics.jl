module spiking
    import Munge
    using StatsBase

    function nonlocality(times::Vector, beh::DataFrame, R::ReceptiveField)
    end
    function nonlocality(beh_of_spikes::Array, R::ReceptiveField)
    end

    function torate(spikes::DataFrame, dims=:unit; kws...)
        torate.(Munge.tensorize(spikes, dims, :time))
    end
    function torate(times::Vector; grid, windowsize)
        # use sample radius from adaptive module
    end
    function torate(times::Vector; grid)
        count = fit(Histogram, grid).weights
        count ./ (grid[2]-grid[1])
    end

end

