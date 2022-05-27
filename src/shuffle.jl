module shuffle

    using Distributions
    using DataFrames
    using DrWatson
    export filt
    include(srcdir("filt.jl"))
    include(srcdir("utils.jl"))
    SplitType = Union{Vector{Symbol}, Symbol, String, Vector{String}}
    defaultFilters = merge(filt.speed_lib)
    using Infiltrator

    # # # # # # # # # # # # # # # # # # # # # # # # 
    # Main shuffle types
    # # # # # # # # # # # # # # # # # # # # # # # # 
    function all(spikes::DataFrame, distribution::Distribution;
            only_times::Bool=false, kws...)::DataFrame
        @debug "all where distribution::Distribution=$distribution"
        #spikes = copy(spikes)
        spikes = transform(spikes, :time => identity => :time, copycols=false)
        spikes.time .+= rand(distribution, 1) 
        only_times ? spikes[!,[:time]] : spikes
    end
    function all(spikes::SubDataFrame, distribution::Distribution; kws...)::SubDataFrame
        @debug "all where distribution::Distribution=$distribution"
        spikes.time .+= rand(distribution, 1) 
        spikes
    end
    function all(spikes::DataFrame; kws...)::DataFrame
        @debug "all where no distribution"
        distribution = _create_distribution(spikes; kws...)
        spikes = all(spikes, distribution)
    end

    function by(spikes::DataFrame, distribution::Distribution; 
                split::SplitType=:unit, sort::Bool=false, kws...)::DataFrame
        @debug "by() with distribution=$distribution"
        spikes = transform(spikes, :time => identity => :time, copycols=false)
        spikes = combine(groupby(spikes, split, sort=sort), sp->all(sp, distribution; kws...))
    end

    function by(spikes::DataFrame; split::SplitType=[:unit,:traj], kws...)::DataFrame
        @debug "by() with no distribution"
        distribution = _create_distribution(spikes; kws...)
        by(spikes, distribution; split=split, kws...)
    end

    function byspike(spikes::DataFrame; kws...)::DataFrame
        distribution = _create_distribution(spikes; kws...)
        byspike(spikes, distribution)
    end

    function byspike(spikes::DataFrame, distribution::Distribution)::DataFrame
        spikes = transform(spikes, :time => identity => :time, copycols=false)
        jitter = rand(distribution, size(spikes, 1))
        spikes.time .+= jitter
        return spikes
    end


    # # # # # # # # # # # # # # # # # # # # # # # # 
    # Internal distribution functions
    # # # # # # # # # # # # # # # # # # # # # # # # 
    function _create_distribution(spikes, distribution::Symbol=:uniform;
            width=:traj, kws...)::Distribution
        if :shuffledist_df in keys(kws)
            @debug ":shuffledist_df in keys"
            data = kws[:shuffledist_df]
        elseif :data in keys(kws)
            data = kws[:data]
        else
            @debug ":shuffledist_df NOT in keys"
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
                               filters::AbstractDict=defaultFilters,
                               kws...)::Real
        G = try
            r(x)    = replace(x, missing=>NaN)
            inds    = (!).(isnan.(r(data.traj))) 
            for (col, filtFunc) in filters
                inds .&= filtFunc(r(data[!, col]))
            end
            data = data[inds, :]
            G = combine(groupby(data, :traj, sort=false),
                              :time=> (x->diff([minimum(x), maximum(x)])) => :range).range;
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
    
    function _session(data::DataFrame)::Real
        session_length = utils.dextrema(data.time)[1]
        @debug "session_length=$session_length"
        return session_length
    end

    # Prepackaged description of shuffles
    standard_shuffles = Dict(
        :cDt_t => (;name="σ(cellᵢ,Δtraj)|trajⱼ",
                    desc=s"""randomly uniform shift a cell on the timescale of
                    a trajectory length"""
                   )
        )

end
export shuffle
