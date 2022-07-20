module metrics

    using StatsBase
    using Entropies: Probabilities
    using NaNStatistics
    using DataFrames
    import ..Field
    import ..Field: ReceptiveField
    import Table
    import Utils
    import Load.utils: register
    import Load: keep_overlapping_times
    using Infiltrator

    export Metrics
    export bitsperspike, bitspersecond, coherence, totalcount, maxrate,
           maxcount 

    mutable struct Metrics
        data::AbstractDict{Symbol, Any}
        Metrics()  = new(Dict{Symbol,Any}())
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
    function Table.to_dataframe(M::Metrics; kws...) 
        kws = (;kws..., key_name="metric")
        Table.to_dataframe(M.data; kws...)
    end

    function push_metric!(R::ReceptiveField, F::Function; 
            name::Union{Symbol,Nothing}=nothing)
        name = name === nothing ? Symbol(F) : name
        push!(R.metrics, name=>F(R))
    end
    pop_metric!(R::ReceptiveField, name::Symbol) = pop!(R.metrics, name)

    skipnan(x) = Iterators.filter(!isnan, x)
    function _convert_to_prob(behProb::AbstractArray)
        if behProb isa AbstractArray
            behProb = Probabilities(collect(skipnan(vec(behProb))))
        end
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


    function bitsperspike(F::T where T<:ReceptiveField)
        bitsperspike(F.rate, F.occ.prob)
    end
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
        beh, spikes = keep_overlapping_times(beh, spikes)
    end

    function totalcount(F::ReceptiveField)
        sum(F.count)
    end
    function maxcount(F::ReceptiveField)
        maximum(F.count)
    end
    function maxrate(F::ReceptiveField)
        nanmaximum(F.rate)
    end
    function meanrate(F::ReceptiveField)
        nanmean(F.rate)
    end
    function convexhull(F::ReceptiveField)
    end

    """
        coherence(F::Field.ReceptiveField)

    computes the spatial coherence of a field
    """
    function coherence(F::ReceptiveField)
        summation = Vector{Float32}(undef, length(vec(F.rate)))
        interaction_mat = get_interaction_mat(ndims(F.rate))
        for (i,ind) in enumerate(eachindex(F.rate))
            subject = F.rate[ind]
            interactions = collect(Tuple(ind)) .+ interaction_mat
            interactions[:,1] = max(min(interactions[:,1], 1), size(F.rate,1))
            interactions[:,2] = max(min(interactions[:,2], 1), size(F.rate,2))
            int_samples = F.rate[eachrow(interactions)]
            summation[i] = sum(int_samples * subject)
        end
    end

    function get_interaction_mat(n::Int)
        poss = [1,-1]
        all_possible_neighbors = Iterators.product((poss for i in 1:n)...)
        all_possible_neighbors = collect(all_possible_neighbors)
        all_possible_neighbors = vcat(all_possible_neighbors...)
        vcat([collect(x)[Utils.na, :] for x in all_possible_neighbors]...)
    end

                                                                
# ,---.               o                            |              
# |---',---.,---..   ..,---.,---.,---.    ,---..  ,|--- ,---.,---.
# |  \ |---'|   ||   |||    |---'`---.    |---' >< |    |    ,---|
# `   ``---'`---|`---'``    `---'`---'    `---''  ``---'`    `---^
#               |                                                 
#                 
#               o     ,---.     
#               .,---.|__. ,---.
#               ||   ||    |   |
#               ``   '`    `---'
#               beyond just the fields

trajdiversity(spikes::AbstractDataFrame, cell::Int) = length(unique(spikes.traj))
function trajdiversity(spikes::DataFrame, beh::DataFrame)
    _, spikes = register(beh,spikes;transfer=["traj"])
    combine(groupby(spikes, :traj), trajdiversity)
end

# CELL COFIRING
function xcorr(rate::Array, cell1::Int, cell2::Int; lags=-200:200)
    x, y = rate[:,1], rate[:,2]
    StatsBase.crosscorr(x,y,lags)
end
function xcorr(spikes::DataFrame)
    units1 = unique(spikes.unit)
    units2 = unique(spikes.unit)
    results = []
    for (cell1,cell2) in Iterators.product(units1,units2)
        push!(results,xcorr(spikes, cell1, cell2))
    end
end


end
