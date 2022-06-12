module recon

    include("../utils.jl")
    include("./operation.jl")  #valid if a module has already nested this!
    using NaNStatistics

    """
    reconstruction(jointdensity, marginalrate) -> R̂

    Example: reconstruct the response to `p(x, y, Γ)` under marginal `p(Γ)`
    """
    function reconstruction(jointdensity::AbstractArray,
            marginalrate::AbstractArray; 
            marginalize_dims::Union{Nothing, Vector{Int}}=nothing, kws...)
        singular_marginal_dims = findall(size(marginalrate).==1)
        if isempty(singular_marginal_dims)
            @warn "Potentially missing singular marginal dims"
        end
        if marginalize_dims == nothing
            dims = setdiff(1:ndims(jointdensity), singular_marginal_dims)
        else
            dims = marginalize_dims
        end
        RateProb = jointdensity .* marginalrate
        Prob = jointdensity
        for dim ∈ dims
            RateProb = nansum(RateProb; dims=dim)
            Prob     = nansum(Prob;     dims=dim)
        end
        return utils.squeeze(RateProb./Prob)
    end

    #"""
    #reconstruction(jointdensity, marginalcount, marginalocc)

    #version for spike count based reoncstruction
    #"""
    #function reconstruction(jointdensity::AbstractArray,
    #        marginalcount::AbstractArray, marginalocc::AbstractArray;
    #        extra_marginal_dims::Vector{Int}=[], kws...)
    #    singular_marginal_dims = findall(size(marginalcount).==1)
    #    singular_marginal_dims = singular_marginal_dims ∪ extra_marginal_dims
    #    if isempty(singular_marginal_dims)
    #        @warn "Potentially missing singular marginal dims"
    #    end
    #    dims = setdiff(1:ndims(jointdensity), singular_marginal_dims)
    #    CountProb = jointdensity .* marginalcount
    #    OccProb   = jointdensity .* marginalocc
    #    for dim ∈ dims
    #        CountProb = nansum(CountProb; dims=dim)
    #        OccProb   = nansum(OccProb;   dims=dim)
    #    end
    #    return utils.squeeze(CountProb./OccProb)
    #end

    """
    reconstruction_error(joint, reconstructed_joint) -> ϵ
    """
    function reconstruction_error(original::Dict, 
                                  reconstructed::Dict;
                                  kws...)
        func(x::AbstractArray, y::AbstractArray) = reconstruction_error(x,y;kws...)
        #println("func->",func)
        operation.apply(func, original, reconstructed)

    end
    function reconstruction_error(field::AbstractArray, 
                                  reconstructed_field::AbstractArray;
                                  L::Int=2, reduce::Bool=true, 
                                  debug::Bool=false)
        ex = extrema(utils.skipnan(vec(field)))
        Δ = diff([ex...])
        ε = vec(field) .- vec(reconstructed_field)
        if all(isnan.(ε))
            @warn "All differnces are NaN"
        end
        if L != 1
            ε  = ε .^ L
        end
        if all(isnan.(ε))
            @warn "(2)All differnces are NaN"
        end
        if reduce
            ε = nanmean(ε)
            ε /= Δ
        else
            ε ./= Δ
        end
        if debug
            p = Plots.plot(
                           Plots.heatmap(utils.squeeze(field), title="marginal"),
                           Plots.heatmap(utils.squeeze(reconstructed_field), title="recon. marginal"),
                title="error=$ε")
            return p
        end
        return ε
    end



end
