quickactivate("/home/ryoung/Projects/goal-code/")
using DataFramesMeta
@time include(scriptsdir("fields", "Include.jl"))
using Plots
using StatsPlots
using LazyGrids: ndgrid
animal, day = "RY16", 36
epoch = 7
task = raw.load_task(animal, day)
boundary = task[(task.name.=="boundary") .& (task.epoch .== epoch), :]
append!(boundary, DataFrame(boundary[1,:]))
_, spikes = raw.register(beh, spikes, on="time", transfer=["velVec","x","y"])

# Get occupancy and fields
include(scriptsdir("fields","Initialize.jl"))
include(scriptsdir("fields","PlaceField.jl"))

# Compute centroids (xy spiking medians)
units = combine(groupby(spikes, :unit), 
                :x=>median, :y=>median,
                :x=>mean,   :y=>mean
               )
function base_plot()
    p = heatmap(place.cgrid..., place.occR', colorbar_title="occupancy", c=:oslo,
                 colorbar_title_location=:top, colorbar_titlefontrotation=0)
    @df beh[1:5:end,:] plot!(:x, :y; alpha=0.2, c=:white, label="pos")
    plot!(boundary.x, boundary.y, c=:orange, linestyle=:dash, label="boundary")
    p
end

p1 = base_plot()
@df units scatter!(:x_median, :y_median, c=:red, label="Cell field medians\n(no vel filter)", ylim=(10,100),
                   alpha=0.5, legend_position=:topleft)


spikes_filt = @subset(spikes, abs.(:velVec) .> 1)
units = combine(groupby(spikes_filt, :unit), 
                :x=>median, :y=>median,
                :x=>mean,   :y=>mean
               )

p2 = base_plot()
@df units scatter!(:x_median, :y_median, c=:red, label="Cell field medians\n(\$v\$ > 1cm/s)", ylim=(10,100),
                   alpha=0.5,
                  legend_position=:topleft, legend_fontpointsize=30)


speed=2
spikes_filt = @subset(spikes, abs.(:velVec) .> speed)
units = combine(groupby(spikes_filt, :unit), 
                :x=>median, :y=>median,
                :x=>mean,   :y=>mean
               )
p3 = base_plot()
@df units scatter!(:x_median, :y_median, c=:red, label="Cell field medians\n(\$v\$ > 1cm/s)", ylim=(10,100),
                   alpha=0.5,
                  legend_position=:topleft, legend_fontpointsize=30)
plot!(boundary.x, boundary.y, c=:black, label="boundary")
@df units scatter!(:x_median, :y_median, c=:red, label="Cell field medians\n(\$v\$ > $(speed)cm/s)", ylim=(10,100),
                   alpha=0.5,
                  legend_position=:topleft, legend_fontpointsize=30)

speed=5
spikes_filt = @subset(spikes, abs.(:velVec) .> speed)
units = combine(groupby(spikes_filt, :unit), 
                :x=>median, :y=>median,
                :x=>mean,   :y=>mean
               )
p4 = base_plot()
@df units scatter!(:x_median, :y_median, c=:red, 
                   label="Cell field medians\n(\$v\$ > $(speed)cm/s)", ylim=(10,100),
                   alpha=0.5,
                  legend_position=:topleft, legend_fontpointsize=30)

p = Plots.plot(p1,p4, size=(1200,600), dpi=400)
utils.savef("cells","coverage_xy_medians")

# Compute centroids (medians per field)
K = [keys(place.Rₕ)...]
X, grid = place.Cₕ[K[1]], place.cgrid

function get_centroid(X; grid=nothing, func=nanmean, upper=0.999, lower=0.001)
    grid = ndgrid(grid...)
    lower = nanquantile(X, lower)
    upper = nanquantile(X, upper)
    X[X .< lower] .= NaN
    X[X .> upper] .= NaN
    [func(vec(X.*g))./func(vec(X)) for g in grid]
end
Y = hcat([get_centroid(place.Cₕ[K[k]]; grid=grid, func=nansum) for k in 1:length(K)]...)'
p5 = base_plot()
scatter!(Y[:,1], Y[:,2], c=:red, label="Cell field medians\n(\$v\$ > 2cm/s)",
         ylim=(10,100), alpha=0.5, legend_position=:topleft,
         legend_fontpointsize=30)
utils.savef("cells","coverage_occupancy_norm_centroids")


# Color convex hulls of fields
using LazySets
function hull_of_field(X, grid; thresh=0.95)
    thresh = nanquantile(vec(X), thresh)
    bool = X .> thresh
    G = ndgrid(grid...)
    G = [vec(g[bool]) for g in G]
    points = [[p...] for p in zip(G...)]
    h = ConvexHullArray([Ball2(p, 0.0) for p in points])
    h, bool
end


