module Shuf

    __revise_mode__ = :evalassign
    using Distributions
    using DataFrames
    using DrWatson

    using Infiltrator
    import Shuffle: shuffle!

    using Filt
    SplitType = Union{Vector{Symbol}, Symbol, String, Vector{String}}
    defaultFilters = merge(Filt.speed_lib)
    using Utils

    distribution_based = [:addSampleOfDist, :jitterBy, :jitterAllSpikes]
    function isaDistributionShuffle(shuffle_func::Symbol)
        any([shuffle_func === func for func in distribution_based])
    end
    function isaDistributionShuffle(shuffle_func::Function)
        isaDistributionShuffle(Symbol(shuffle_func))
    end

    # ===================================
    # Prepackaged description of shuffles
    # ===================================

    """
    standard_shuffles

    link shortcut namess of shuffles to a description and set of arguments
    controlling the application of shuffles with thee machinery below

    | Input | Description |
    |:------|:------------|
    |fullname     | A descritive fullname that you would call it, space-allowing|
    |desc         | either a textual or mathematical description|
    |shuffle_func | the top-level function called to shuffle|
    |prop         | what property/column to jitter or permute|
    |split        | (optional) split jittering or permuting by this?|
    |distribution | (optional) if jitter, what distribution do we draw from?|

    """
    standard_shuffles = Dict(

        :cDt_t => (;fullname="+=x~σ(cellᵢ,Δtraj)|trajⱼ",
                    desc=s"""randomly uniform shift a cell on the timescale of
                    a trajectory length""",
                    shuffle_func = :jitterBy,
                    split        = [:unit, :traj],
                    distribution = :uniform,
                    prop         = :time
                   ),

        :dotson => (;fullname="split by trajreltime_bin and permute traj",
                    desc=s"""for each trajectory, set the time on a scale 0 to
                    1, where 0 and 1 are beginning and end. shuffle spikes to
                    keep the their value i ∈ [0,1].""",
                    shuffle_func = :permuteBy,
                    split        = [:unit, :trajreltime_bin],
                    prop         = :time
                   ),
       )

    """
    `applyStandardShuffle`

    Tool for applying standard shuffle *presets*. You can make
    custom shuffles with the jitty and permute functions below
    if you have tidy data frames.

    # Inputs

    | Inputs | Description |
    |:-------|:------------|
    | name   | symbolic name of a [`standard_shuffles`](@ref) preset|

    # Outputs
    | Outputs | Description |
    |:--------|:------------|
    |partial_dist

    # Downstream pos/kws set by presets
    ## jitterBy
    split, prop
    ### ```_create_distribution```
    distribution, width, shuffledist_df, data
    #### ```_```**width```_```function**
    variable sets of arguments
    ## permuteBy
    split, prop
    """
    function applyStandardShuffle(name::Symbol; precompute_dist::Bool=true)
        settings = standard_shuffles[name]
        applyStandardShuffle(settings; precompute_dist)
    end
    function applyStandardShuffle(settings::AbstractDict; precompute_dist::Bool=true)
        settings = utils.namedtuple_to_dict(settings)
        applyStandardShuffle(settings; precompute_dist)
    end
    function applyStandardShuffle(settings::NamedTuple; precompute_dist::Bool=true)
        @info "Settings=$settings"
        dist_type = occursin("jitter", lowercase(String(settings[:shuffle_func])))
        func = eval(settings[:shuffle_func])
        if dist_type && precompute_dist
            partial_dist(spikes::DataFrame) = begin
                @info "distribution_settings=$settings"
                _create_distribution(spikes::DataFrame, settings[:distribution];
                                     settings...)
            end
        end
        partial(pos...; newkws...) = func(pos...; settings..., newkws...)
        if dist_type
            return partial_dist, partial
        else
            return partial
        end
    end

    # # # # # # # # # # # # # # # # # # # # # # # # 
    # Atomic shuffle operations
    # # # # # # # # # # # # # # # # # # # # # # # # 
    

    # ---- Add a sampled value from a distribution to a property ----- 
    function addSampleOfDist(spikes::DataFrame, distribution::Distribution;
            only_prop::Bool=false, kws...)::DataFrame
        @debug "addSampleOfDist where distribution::Distribution=$distribution"
        #spikes = copy(spikes)
        spikes = transform(spikes, :time => identity => :time, copycols=false)
        spikes.time .+= rand(distribution, 1) 
        only_prop ? spikes[!,:time] : spikes
    end
    function addSampleOfDist(spikes::SubDataFrame, distribution::Distribution; 
            kws...)::SubDataFrame
        spikes.time .+= rand(distribution, 1) 
        spikes
    end
    function addSampleOfDist(spikes::DataFrame; kws...)::DataFrame
        @debug "addSampleOfDist where no distribution"
        distribution = _create_distribution(spikes; kws...)
        spikes = addSampleOfDist(spikes, distribution)
    end
    
    # ---- Permute samples of a property -----
    function permuteSamples(spikes::DataFrame; prop::Symbol=:time, kws...)
        spikes = transform(spikes, :time => identity => :time, copycols=false)
        shuffle!(spikes[!, prop])
        spikes
    end

    function permuteSamples(spikes::SubDataFrame; prop::Symbol=:time, 
            keepold::Bool=false, kws...)
        if keepold
            spikes[!,String(prop)*"old"] = spikes[!,prop]
        end
        shuffle!(spikes[!, prop])
        spikes
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
        jitterAllSpikes(spikes, distribution)
    end


    # # # # # # # # # # # # # # # # # # # # # # # # 
    # How to apply shuffles over groups
    # # # # # # # # # # # # # # # # # # # # # # # # 

    # By methods shuffle within a group, preserving that grouping property (could within traj, or a bin# of relative trial time, etc)
    
    """
    DEFINED distribution, split by something 
    """
    function jitterBy(spikes::DataFrame, distribution::Distribution; 
            split::SplitType=[:unit,:traj], sort::Bool=false, kws...)::DataFrame
        @debug "by() with distribution=$distribution"
        @info "jitterBy split=$split"
        spikes = transform(spikes, 
                           :time => identity => :time, 
                           copycols=false)
        spikes = combine(groupby(spikes, split, sort=sort),
                         sp->addSampleOfDist(sp, distribution; kws...))
    end

    """
    UNDEFINED dist, split by something
    """
    function jitterBy(spikes::DataFrame;
            split::SplitType=[:unit,:traj], kws...)::DataFrame
        @debug "by() with no distribution"
        distribution = _create_distribution(spikes; kws...)
        jitterBy(spikes, distribution; split=split, kws...)
    end

    """
    UNDEFINED dist, split by something
    """
    function permuteBy(spikes::DataFrame; sort::Bool=false,
            split::SplitType=[:unit,:traj], kws...)::DataFrame
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
export Shuf
