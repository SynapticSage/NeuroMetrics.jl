
# -----
# Notes
# -----
#
# - At the momemtn I run FullVideoDecode before this script.
# - Baise measures bias AWAY from
#
# MIGHT BE SOME THOUGHT :mind: BUGS HERE
#
# draw the math out


using StatsBase
@assert issorted(cells.unit)
lfp, spikes = raw.register(lfp, spikes, on="time", transfer="phase",
                           tolerance=1/1000, tolerance_violation=NaN)
beh, spikes = raw.register(beh, spikes, on="time", transfer="trajreltime",
                           tolerance_violation=NaN, tolerance=0.04)
cells, spikes = raw.register(cells, spikes, on="unit", transfer=colorby)
rename!(spikes, colorby=>:tau)

# Pyr
cells_pyr = @subset(cells, :meanrate .< 6)
spikes_pyr = spikes[cells.meanrate[spikes.unit] .< 6, :]
cells_pyr, spikes_pyr = raw.cell_resort(cells_pyr, spikes_pyr, :unit)

# =============================================
phase,trajreltime,tau = (spikes_pyr.phase,  spikes_pyr.trajreltime, 
                         cells_pyr[spikes_pyr.unit, colorby])

reltime_vs_tau = combine(groupby(spikes_pyr, :tau), 
                         :trajreltime => median,
                         :trajreltime => mean)
@df reltime_vs_tau Plots.scatter(:trajreltime_mean, :tau)

pmed = @df reltime_vs_tau Plots.scatter(:trajreltime_median, :tau, xlim=(0,1))
pdist = Plots.density(beh.trajreltime, xlabel="Relative traj", xlim=(0,1))
Plots.plot(pmed, pdist, layout=Plots.grid(2,1))

# =============================================
spikes_pyr.binned_reltime = utils.searchsortednearest.([0:0.05:1], spikes_pyr.trajreltime)

tau_vs_reltime = combine(groupby(spikes_pyr, :binned_reltime), 
                         :tau => x->mean(utils.skipnan(x)),
                         :tau => nanmedian)
@df tau_vs_reltime Plots.plot(:binned_reltime, :tau_function, alpha=0.6, label="mean tau @ relative traj bin")
Plots.hline!([0], c=:black, linestyle=:dash, label=nothing)
Plots.hline!([nanmaximum(spikes_pyr.tau)], c=:red, linestyle=:dash,label=nothing)
Plots.hline!([nanminimum(spikes_pyr.tau)], c=:red, linestyle=:dash,label=nothing)

pmed = @df reltime_vs_tau Plots.scatter(:trajreltime_median, :tau, xlim=(0,1))
pdist = Plots.density(beh.trajreltime, xlabel="Relative traj", xlim=(0,1))
Plots.plot(pmed, pdist, layout=Plots.grid(2,1))

V = []
groups = groupby(spikes_pyr, :binned_reltime)
for group in groups
    group = dropmissing(group)
    x = collect(utils.skipnan(group.tau))
    push!(V,Plots.violin(x))
end
