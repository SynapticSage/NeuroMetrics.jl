module fit

using LsqFit
using Statistics

_vonmises(γᵢ, C₁,C₂, κ, γ) = C₁ .* exp.(κ.*cos.(γ .- γᵢ)) .+ C₂
_vonmises(x,p) = _vonmises(x, p[1],p[2],p[3],p[4])

_gaussian(xᵢ, κ, μₓ, σₓ, C₁) = κ * pdf(Gaussian(μₓ, σₓ), xᵢ) .+ C₁
_gaussian(x,p) = _gaussian(x, p[1],p[2],p[3],p[4])

function fitvonmises(field, grid)
    fit = fit_curve(_vonmises, grid, field)
    ideal = _vonmises(grid, fit.params...)
    (fit=fit, params=fit.params, ideal=ideal)
end

function fitgaussian(field, grid)
    fit = fit_curve(_gaussian, grid, field)
    ideal = _gaussian(grid, fit.params...)
    (fit=fit, params=fit.params, ideal=ideal)
end
    
module metric
    function rayeligh(field, grid)
    end
    function maxmetric(field, grid)
    end
end

end
