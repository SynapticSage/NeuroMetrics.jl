import Plots
includet(srcdir("utils.jl"))
import .utils

# Load cycles with important goodies


# WANT LIST
# - Theta direction vectors
#   - color points by relative traj time
#   - color/split points by 
#       - stopwell
#       - futurestopwell
#       - paststopwell
# - 


colorby = splitby = nothing

vector = collect(utils.skipnan(cycles.curfinal))
vector = vector[abs.(vector) .< quantile(abs.(vector), [0.99999])]

# ANGULAR BINS
θ=-π:0.1:π; 
θᵪ = mean([θ[begin:end-1] θ[begin+1:end]];dims=2); 
r=fit(Histogram, vector, θ).weights; 

hist_1 = Plots.plot(θᵪ, r, proj=:polar, 
           label="ThetaSequence direction", 
           title="Most of the theta sequences point away from home", 
           ylim=(0, 175), fill=0,
           alpha=0.5)
utils.savef("theta","goal","hist_sequences_to_and_away_from_home")

hist_2 = Plots.stephist(angle.(vector), proj=:polar, 
           label="ThetaSequence direction", 
           title="Most of the theta sequences point away from home", 
           ylim=(0, 175), fill=0,
           alpha=0.5)
utils.savef("theta","goal","stephist_sequences_to_and_away_from_home")

hist_combo = Plots.plot(hist_1, hist_2, title="")

# RAW VECTORS IN SCATTER FORM
pol_1 = Plots.scatter(angle.(vector), abs.(vector), proj=:polar,
              ylim=quantile(abs.(vector), [0, 0.99]),
           label="ThetaSequence direction+length", 
           title="", 
           markerstrokewidth=0,
           markersize=3,
           alpha=0.35)
pol_2 = Plots.scatter(angle.(vector), abs.(vector), proj=:polar,
              ylim=quantile(abs.(vector), [0, 0.7]),
           label="ThetaSequence direction+length", 
           title="", 
           markerstrokewidth=0,
           markersize=3,
           alpha=0.25)
pol_3 = Plots.scatter(vector, 
           label="ThetaSequence direction+length", 
           title="", 
           markerstrokewidth=0,
           markersize=3,
           alpha=0.35)
pol_4 = Plots.scatter(vector, 
                      ylim=quantile(abs.(vector), [0.7, 0.7]).*[-1,1],
                      xlim=quantile(abs.(vector), [0.7, 0.7]).*[-1,1],
           label="ThetaSequence direction+length", 
           title="", 
           markerstrokewidth=0,
           markersize=3,
           alpha=0.25)
Plots.plot(pol_1, pol_2, pol_3, pol_4, linewidth=0)

utils.savef("theta","goal",
            "scatter_sequences_to_and_away_from_home-more away from home--but longer to home")


# Statistics of length in each bin
θ = -π : π / 4 : π
D=Dict()
D[:θₙ] = utils.searchsortednearest.([θ], angle.(vector))
D[:r]  = abs.(vector)
D = DataFrame(D)
D = combine(groupby(D, :θₙ), :r=>mean, :r=>sum, :r=>median)
D[!,:θ] = θ[D[!,:θₙ]]
thetar_mean   = @df D Plots.plot(:θ, :r_mean,   label="mean",   proj=:polar, fill=0, alpha=0.5)
thetar_median = @df D Plots.plot(:θ, :r_median, label="median", proj=:polar, fill=0, alpha=0.5)
thetar_sum    = @df D Plots.plot(:θ, :r_sum,    label="sum",    proj=:polar, fill=0, alpha=0.5)
Plots.plot(thetar_mean, thetar_median, thetar_sum)
utils.savef("theta","goal",
            "stats_of_theta_vectorlengths_per_bin")


