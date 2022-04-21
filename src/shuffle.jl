module shuffle

    using Distributions
    using DataFrames
    using DrWatson
    export filt
    include(srcdir("filt.jl"))
    include(srcdir("utils.jl"))
    SplitType = Union{Vector{Symbol}, Symbol, String, Vector{String}}
    defaultFilters = merge(filt.speed_lib)

    # # # # # # # # # # # # # # # # # # # # # # # # 
    # Main shuffle types
    # # # # # # # # # # # # # # # # # # # # # # # # 
    function all(spikes::AbstractDataFrame, distribution::Distribution;
            only_times::Bool=false, kws...)
        @debug "all where distribution::Distribution=$distribution"
        spikes = copy(spikes)
        spikes.time = spikes.time .+ rand(distribution, 1) 
        only_times ? spikes[!,[:time]] : spikes
    end
    function all(spikes; kws...)
        @debug "all where no distribution"
        spikes = copy(spikes)
        distribution = _create_distribution(spikes; kws...)
        spikes = all(spikes, distribution)
    end

    function by(spikes::DataFrame, distribution::Distribution; 
                split::SplitType=:unit, kws...)
        @debug "by() with distribution=$distribution"
        spikes = copy(spikes)
        spikes = combine(groupby(spikes, split), sp->all(sp, distribution; kws...))
    end

    function by(spikes; split::SplitType=:unit, kws...)
        @debug "by() with no distribution"
        distribution = _create_distribution(spikes; kws...)
        by(spikes, distribution; split=split, kws...)
    end

    function byspike(spikes, distribution::Distribution)
        spikes = copy(spikes)
        nSpikes = length(spikes)
        if distribution isa String
        end
        jitter = rand(distribution, nSpikes)
        spikes.time .+= jitter
        return spikes
    end

    function byspikes(spikes; kws...) 
        distribution = _create_distribution(spikes; kws...)
        byspike(spike, distribution)
    end

    # # # # # # # # # # # # # # # # # # # # # # # # 
    # Internal distribution functions
    # # # # # # # # # # # # # # # # # # # # # # # # 
    function _create_distribution(spikes, distribution::Symbol=:uniform; width=:traj, kws...)
        if :shuffledist_df in keys(kws)
            @debug "data in keys"
            data = kws[:shuffledist_df]
        else
            @debug "data NOT in keys"
            data = spikes
        end
        if width == :traj
            width = _typical_trajtime(data; kws...)
        elseif width == :session
            width = _session(data)
        end
        if distribution == :uniform
            return Uniform(-width/2, width/2)
        elseif distribution == :normal
            return Normal(0, width)
        else
            throw(ArgumentError("pass a distribution or an accepted string for distribution. curren val=$distribution"))
        end
    end

    # # # # # # # # # # # # # # # # # # # # # # # # 
    # WIDTH Functions
    # # # # # # # # # # # # # # # # # # # # # # # # 
    function _typical_trajtime(data::DataFrame;
                               filters::AbstractDict=defaultFilters, kws...)
        G = try
            inds    = (!).(isnan.(data.traj)) 
            for (col, filtFunc) in filters
                inds .&= filtFunc(data[!, col])
            end
            data = data[inds, :]
            G = combine(groupby(data, :traj),
                      :time=> (x->diff([extrema(x)...])) => :range).range;
        catch e
            if e isa ArgumentError
                @warn "ArugmentError: try using data=beh if keys are missing"
            end
            throw(e)
        end
        median_traj_size = median(G)
        @debug "median_traj_size=$median_traj_size"
        return median_traj_size
    end
    
    function _session(data::DataFrame)
        session_length = utils.dextrema(data.time)
        @debug "session_length=$session_length"
        return session_length
    end

end
export shuffle
