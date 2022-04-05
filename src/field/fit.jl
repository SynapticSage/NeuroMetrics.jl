module fit

using LsqFit
using Statistics

_vonmises(γᵢ, C₁,C₂, κ, γ) = C₁ .* exp.(κ.*cos.(γ .- γᵢ)) .+ C₂
_vonmises(x,p) = _vonmises(x, p[1],p[2],p[3],p[4])

_gaussian(xᵢ, κ, μₓ, σₓ, C₁) = κ * pdf(Gaussian(μₓ, σₓ), xᵢ) .+ C₁
_gaussian(x,p) = _gaussian(x, p[1],p[2],p[3],p[4])

function fitvonmises(field::AbstractArray, grid::AbstractArray)
    fit = fit_curve(_vonmises, grid, field)
    ideal = _vonmises(grid, fit.params...)
    (fit=fit, params=fit.params, ideal=ideal)
end

function fitgaussian(field::AbstractArray, grid::AbstractArray)
    fit = fit_curve(_gaussian, grid, field)
    ideal = _gaussian(grid, fit.params...)
    (fit=fit, params=fit.params, ideal=ideal)
end
    
module metric
    """
    rayleigh

    field could be the empirical field or the `fitted field tuning curve`

    recommend the smoothed/fitted
    """
    function rayeligh(field::AbstractArray, grid::AbstractArray)
        n = length(grid)
        π/(n*sin(π/n)) * (sum(field.*grid)/sum(field))
    end
    """
    maxmetric

    used for PDI and EDI pg15, sarel et al supplemental

    field could be the empirical field or the `fitted field tuning curve`

    recommend the smoothed/fitted
    """
    function maxmetric(field::AbstractArray)
        maximum(field)/mean(field)
    end
end

end
