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
    import Shuffle: shuffle!

    # ===================================
    # Prepackaged description of shuffles
    # ===================================
    standard_shuffles = Dict(
        :cDt_t => (;name="spiketime += x ~ σ(cellᵢ,Δtraj)|trajⱼ",
                    desc=s"""randomly uniform shift a cell on the timescale of
                    a trajectory length""",
                    by=nothing, 
                    distribution=nothing,  
                    prop=:time)
       )
    # by: a grouping that the distribution applied within
    # distribution: a distribution we sampled from an added to the prop or permutation of the prop,
    # prop: the property that is jittered or permuted)

    #function applyStandardShuffle()
    #end

    # # # # # # # # # # # # # # # # # # # # # # # # 
    # Atomic shuffle operations
    # # # # # # # # # # # # # # # # # # # # # # # # 
    

    # ---- Add a sampled value from a distribution to a property ----- 

    function addSampleOfDist(spikes::DataFrame, distribution::Distribution;
            only_prop::Bool=false, prop::Symbol=:time, kws...)::DataFrame
        @debug "addSampleOfDist where distribution::Distribution=$distribution"
        #spikes = copy(spikes)
        spikes = transform(spikes, prop => identity => prop, copycols=false)
        spikes.time .+= rand(distribution, 1) 
        only_prop ? spikes[!,[prop]] : spikes
    end
    function addSampleOfDist(spikes::SubDataFrame, distribution::Distribution; 
            prop::Symbol=:time, kws...)::SubDataFrame
        @debug "addSampleOfDist where distribution::Distribution=$distribution"
        spikes[!,prop] .+= rand(distribution, 1) 
        spikes
    end
    function addSampleOfDist(spikes::DataFrame; kws...)::DataFrame
        @debug "addSampleOfDist where no distribution"
        distribution = _create_distribution(spikes; kws...)
        spikes = addSampleOfDist(spikes, distribution)
    end
    
    # ---- Permute samples of a property -----

    function permuteSamples(spikes::DataFrame; prop::Symbol=:time, kws...)
        spikes = transform(spikes, prop => identity => prop, copycols=false)
        shuffle!(spikes[!, prop])
    end

    function permuteSamples(spikes::SubDataFrame; prop::Symbol=:time, kws...)
        shuffle!(spikes[!, prop])
    end

    # ----- Jitter per spike elelemnt -----

    """
    non-split distributino readout

    how is this different from all above
    """
    function jitterAllSpikes(spikes::DataFrame, distribution::Distribution)::DataFrame
        spikes = transform(spikes, :time => identity => :time, copycols=false)
        jitter = rand(distribution, size(spikes, 1))
        spikes.time .+= jitter
        return spikes
    end
    """
    undefined dist, no split
    """
    function jitterAllSpikes(spikes::DataFrame; kws...)::DataFrame
        distribution = _create_distribution(spikes; kws...)
        byspike(spikes, distribution)
    end


    # # # # # # # # # # # # # # # # # # # # # # # # 
    # How to apply shuffles over groups
    # # # # # # # # # # # # # # # # # # # # # # # # 

    # By methods shuffle within a group, preserving that grouping property (could within traj, or a bin# of relative trial time, etc)
    
    """
    DEFINED distribution, split by something 
    """
    function jitterBy(spikes::DataFrame, distribution::Distribution; 
                split::SplitType=:unit, sort::Bool=false, kws...)::DataFrame
        @debug "by() with distribution=$distribution"
        spikes = transform(spikes, :time => identity => :time, copycols=false)
        spikes = combine(groupby(spikes, split, sort=sort),
                         sp->addSampleOfDist(sp, distribution; kws...))
    end

    """
    UNDEFINED dist, split by something
    """
    function jitterBy(spikes::DataFrame; split::SplitType=[:unit,:traj], kws...)::DataFrame
        @debug "by() with no distribution"
        distribution = _create_distribution(spikes; kws...)
        by(spikes, distribution; split=split, kws...)
    end
    """
    UNDEFINED dist, split by something
    """
    function permuteBy(spikes::DataFrame; split::SplitType=[:unit,:traj], kws...)::DataFrame
        @debug "by() with distribution=$distribution"
        spikes = transform(spikes, :time => identity => :time, copycols=false)
        spikes = combine(groupby(spikes, split, sort=sort),
                         sp->permuteSamples(sp; kws...))
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
        elseif width == :extrema
            width = _extrema(data; kws...)
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

    function _extrema(data::DataFrame; col=:time)
        extrema(data[!,col])
    end

end
export shuffle
