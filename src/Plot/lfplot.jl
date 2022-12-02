module lfplot

using RecipesBase
using StatsPlots
using DataFrames

"""
    cycle_plot(lfp::DataFrame)

Checks that your lfp dataframe cycle labels are correct (what you expect)
"""
function cycleplot(lf::DataFrame)
    @df lf[1:2500,:] begin
        Plots.plot(:time, :raw, label="raw")
        Plots.plot!(:time, mod2pi.(:phase) .+100,label="phase")
        Plots.plot!(:time, 10*:cycle)
    end
end

@userplot PlotPhaseLock
"""
    plotphaselock(spikes::DataFrame)

Carries out phase locking plot
"""
function plotphaselock(plt::PlotPhaseLock; time=:time)
    spikes = plt.args[1]
    spikes isa DataFrame ? TypeError(:spikes, "1st arg should be dataframe") : nothing
    seriestype --> :histogram
    (spikes[!,time])
end

end
