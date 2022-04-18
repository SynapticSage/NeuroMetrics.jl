module shuffle

    using Distributions

    # # # # # # # # # # # # # # # # # # # # # # # # 
    # Main shuffle types
    # # # # # # # # # # # # # # # # # # # # # # # # 
    function byneuron(spikes)
        nNeurons = unique(spikes.unit)
    end
    function byspike(spikes; distribution::Distribution=Normal())
        spikes = copy(spikes)
        nSpikes = length(spikes)
        if distribution isa String
        end
        jitter = rand(distribution, nSpikes)
        spikes.time .+= jitter
        return spikes
    end
    function byspikes(spikes; distribution::String; kws...)
        distribution = _create_distribution(spikes, distribution; kws...)
    end

    # # # # # # # # # # # # # # # # # # # # # # # # 
    # Internal distribution functions
    # # # # # # # # # # # # # # # # # # # # # # # # 
    function _create_distribution(spikes, distribution::String; width=:traj, kws...)
        if width == :traj
            width = _typical_trajtime(spikes, distribution)
        elseif width == :session
            width = _session(spikes, distribution)
        end
        else
        if distribution == "uniform"
            return Uniform(width)
        elseif distribution == "normal"
            return Normal(0, width)
        else
            throw(ArgumentError("pass a distribution or an accepted string for distribution"))
        end
    end

    # # # # # # # # # # # # # # # # # # # # # # # # 
    # WIDTH Functions
    # # # # # # # # # # # # # # # # # # # # # # # # 
    function _typical_trajtime(data::DataFrame;
                               filters::AbstractDict=merge(filt.speed_lib))
        inds    = (!).(isnan.(beh.traj)) 
        for (col, filtFunc) in filters
            inds .&= filtFunc(data[!, col])
        end
    end
    function _session()
    end

end
export shuffle
