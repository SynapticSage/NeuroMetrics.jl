
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
lfp, spikes = raw.register(lfp, spikes, on="time", transfer="phase")
lfp, spikes = raw.register(lfp, spikes, on="time", transfer="phase")

function plot_phase_by_shift(spikes, cells; shiftcol, plot_title)

    data = hcat(spikes.phase, cells[spikes.unit, shiftcol])
    inds = findall(vec(all((!).(isnan.(data)),dims=2)))
    data = data[inds, :]

    X = fit(Histogram, Tuple((d for d in eachcol(data))))
    phase, shift = [field.edge_to_center(x) for x in X.edges]
    # phase x shift
    X = Float64.(X.weights)
    Xnorm_phase = X ./ nanmaximum(X, dims=1)
    Xnorm_shift = X ./ nanmaximum(X, dims=2)
    Xnorm_phaseshift = X ./ nanmaximum(X, dims=2)./ nanmaximum(X, dims=1)

    clim(item) = Tuple([nanquantile(vec(item), q) for q in (0.27, 0.99)])

   p2= Plots.heatmap(phase, shift, Xnorm_phase',  title="Normalized \nwithin phase", clim=clim(Xnorm_phase))
   Plots.plot!(phase, nanmean(Xnorm_phase .* shift')./nanmean(Xnorm_phase; dims=2), c=:blue, linewidth=3, label="Bias")
   Plots.vline!([0], c=:white, linestyle=:dash, label=nothing)

   p4= Plots.heatmap(phase, shift, Xnorm_shift',  title="Normalized \nwithin shift", clim=clim(Xnorm_shift))
   Plots.plot!(phase, nanmean(Xnorm_shift .* shift')./nanmean(Xnorm_phase; dims=2), c=:blue, linewidth=3, label="Bias")
   Plots.vline!([0], c=:white, linestyle=:dash, label=nothing)

   p3= Plots.heatmap(phase, shift, Xnorm_phaseshift',  title="Normalized \nwithin phase and shift", clim=clim(Xnorm_phaseshift))
   Plots.plot!(phase, nanmean(Xnorm_phaseshift .* shift')./nanmean(Xnorm_phaseshift; dims=2), c=:blue, linewidth=3, label="Bias")
   Plots.vline!([0], c=:white, linestyle=:dash, label=nothing)

   T = Plots.plot(;plot_title, grid = false, showaxis = false, yticks=[],
                  xticks=[],
                 bottom_margin = -50Plots.px) 

    layout = Plots.@layout [A{0.2h}; a b; c _]
    Plots.plot(T, p2, p3,
               Plots.heatmap(phase, shift, X', title="Pure \nspike counts", clim=clim(X)),
               layout=layout
    )
end

# All cells
pall = plot_phase_by_shift(spikes, cells, shiftcol=colorby, plot_title="All cells")



# Interneuron
cell_int = @subset(cells, :meanrate .>= 6)
spikes_int = spikes[cells.meanrate[spikes.unit] .>= 6, :]
cell_int, spikes_int = raw.cell_resort(cell_int, spikes_int, :unit)
pint = plot_phase_by_shift(spikes_int, cell_int, shiftcol=colorby, plot_title="Int")

#Plots.plot(
#   Plots.heatmap(X, title="No axes"),
#   Plots.heatmap(phase, shift, X', title="transpose with axes"),
#   Plots.heatmap(phase, shift, X, title="no transpose with axes"),
#)

