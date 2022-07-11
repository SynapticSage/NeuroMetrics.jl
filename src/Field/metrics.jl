module metrics
    using Entropies: Probabilities
    using NaNStatistics
    import ..Field
    import ..Field: ReceptiveField
    using Infiltrator

    mutable struct Metrics
        data::AbstractDict{Symbol, Any}
        Metrics() = new(Dict{Symbol,Any}())
        Metrics(x) = new(x)
    end
    function Base.getindex(M::Metrics, index...)
        Base.getindex(M.data, index...)
    end
    function Base.setindex!(M::Metrics, val, index::Symbol)
        M.data[index] = val
    end
    Base.push!(M::Metrics, p::Pair{Symbol, <:Any}) = push!(M.data, p)
    Base.pop!(M::Metrics, key)::Any = pop!(M.data, key)
    Base.iterate(M::Metrics) = iterate(M.data)
    Base.iterate(M::Metrics, i::Int64) = iterate(M.data, i)
    Base.length(M::Metrics)  = length(M.data)
    function Base.string(S::T where T<:Metrics; sigdigits=2)
        M = ["$k=$(round(v;sigdigits))" for (k,v) in S]
        join(M, " ")
    end

    skipnan(x) = Iterators.filter(!isnan, x)
    function _convert_to_prob(behProb::AbstractArray)
        if behProb isa AbstractArray
            behProb = Probabilities(collect(skipnan(vec(behProb))))
        end
    end

    """
        coherence(F::Field.ReceptiveField)

    computes the spatial coherence of a field
    """
    function coherence(F::ReceptiveField)
    end
    information(F::T where T <: ReceptiveField; method=:spatialinformation) = information(F.rate, F.occ.prob; method)
    function information(F::Dict, behProb::Probabilities; method=:spatialinformation)
        I = Dict()
        method = eval(method)
        for i in keys(F)
            I[i] = method(F[i], behProb)
        end
        return I
    end
    information(F::Array, behProb::Probabilities; method=:spatialinformation) = eval(method)(F, behProb)


    function bitsperspike(firingrate::AbstractArray, behProb::Probabilities)
        firingrate = vec(firingrate)
        FRoverMeanFR = firingrate ./ nanmean(firingrate)
        # paper r = ∑ pᵢ * rᵢ
        # what I've written here: r = ∑ R_occᵢ / N = μ(R_occᵢ), not the actual rate, but occ norm rate
        I = nansum( behProb .* FRoverMeanFR .* log2.(FRoverMeanFR) )
        # what i've written here: ∑ p(state) * r(state)/r̅ * log₂( r(state) / r̅ )
        # paper : ∑ pᵢ * ( rᵢ/r̅ ) * log₂( rᵢ / r̅ )
        return I
    end
    spatialinformation = bitsperspike
    function bitspersecond(firingrate::AbstractArray, behProb::Probabilities)
        firingrate = vec(firingrate)
        FRoverMeanFR = firingrate ./ nanmean(firingrate)
        # paper r = ∑ pᵢ * rᵢ
        # what I've written here: r = ∑ R_occᵢ / N = μ(R_occᵢ), not the actual rate, but occ norm rate
        I = nansum( behProb .* firingrate .* log2.(FRoverMeanFR) )
        # what i've written here: ∑ p(state) * r(state)/r̅ * log₂( r(state) / r̅ )
        # paper : ∑ pᵢ * ( rᵢ/r̅ ) * log₂( rᵢ / r̅ )
        return I
    end

    function kraskov_mutualinformation(spikeRate, behRate; kws...)
        K = Entropies.Kraskov(;kws...) 
        D = DataSet(spikeRate, behCount)
        HXY = Entropies.genentropy(D, K)
        HX  = Entropies.genentropy(spikeRate, K)
        HY  = Entropies.genentropy(behRate, K)
        HX + HY - HXY
    end
    
    """
    mutualinformation

    get the mutual information of a single cell
    """
    function mutualinformation(beh, spikes, props, unit)
        
        beh, spikes = Load.keep_overlapping_times(beh, spikes)
    end

end
export info
