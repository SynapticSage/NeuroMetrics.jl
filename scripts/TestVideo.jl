# For testing strategies to make videos

using GLMakie
import Plots
includet(srcdir("utils.jl"))


nTime = Int(1e5)
nNeuron = 20
nSpike = Int(3*(nTime/20))
nX, nY = 20, 20

# Create spiking scatter
#sort(collect(eachrow([Float32.(rand(1:nNeuron, spikes)) rand(spikes)])), by=x->x[2])
spikes = hcat(sort(collect(eachrow([Float32.(rand(1:nNeuron, nSpike)) rand(nSpike)])), by=x->x[2])...)'
probability = rand(nX,nY,nTime)
probability ./= sum(probability, dims=(1,2))
t_prob  = LinRange(0,1,nTime)


function select_spikes(t, spikes, time_col=2)
    inds = (spikes[:,time_col] .> t - 0.01) .&& (spikes[:,time_col] .< t+0.01)
    spikes = spikes[inds, :]
    spikes[:,time_col] = spikes[:,time_col] .- t
    spikes
end
function select_prob(t)
    i = utils.searchsortednearest(t_prob, t)
    utils.squeeze(probability[:,:,i])
end

t = Observable(Float64(0))
subset_spikes = @lift select_spikes($t, spikes)
xlim_range    = @lift ($t-0.01, $t+0.01)
scatter_data  = @lift( collect(Point2f(row[end:-1:begin]) for row in eachrow($subset_spikes)) )
fig, ax, sc= scatter(scatter_data, depth_shift=0f0)
xr = LinRange(-0.01, 0.01, size(probability,1))
yr = LinRange(1, 20, size(probability,2))
hm = heatmap!(xr, yr, select_prob(t[]), transparency=true, colormap=(:bamako, 0.2))

on(t) do t
    println("t=$t")
end


framerate = 30
timestamps = range(0, 1, step=0.001)
record(fig, "test.mp4", timestamps) do stamp
    t[] = stamp
    hm[3][] = select_prob(t[])
end
