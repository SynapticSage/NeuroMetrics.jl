"""
    KEY TIMESHIFT FIGURES

Anything that I consider key -- I will try to eventually
migrate here

# Requirements
Run and cache all of the key dictionaries
- I :: Dict of actual measurements per conditions
- S :: Dict of shuffle measurements per conditions
- F :: Dict of field measurements per conditions

# FIGURES

"""


Plot.setfolder(parent,"fixed_heatmaps2")
for key in keys(F)
    f = prep(matrixform(F[key]))
    #push_metric!(f, metrics.bitsperspike)
    #push_shiftmetric!(f, best_tau!; metric=:bitsperspike)
    inds = sortperm(f[:bestshift_bitsperspike][:,1])
    bps  = f[:bitsperspike][inds, :]
    XX = hcat([Utils.norm_extrema(b) for b in eachrow(bps)]...)'
    heatmap(shifts, collect(1:size(XX,1)), XX, title="$(key.datacut)", size=(600,1200))
    vline!([0],c=:black,linestyle=:dash, linewidth=2)
    Plot.save((;xlim="full", key...))
    xlims!(-2,2)
    Plot.save((;xlim=(-2,2), key...))
    xlims!(-1,1)
    Plot.save((;xlim=(-1,1), key...))
end

Plot.setfolder(parent,"fixed_meanmaps2")
P=[]
for key in keys(F)
    f = prep(matrixform(F[key]))
    #push_metric!(f, metrics.bitsperspike)
    #push_shiftmetric!(f, best_tau!; metric=:bitsperspike)
    inds = sortperm(f[:bestshift_bitsperspike][:,1])
    bps  = f[:bitsperspike][inds, :]
    XX = hcat([Utils.norm_extrema(b) for b in eachrow(bps)]...)'
    p=plot(shifts, vec(mean(XX,dims=1)), title="$(key.datacut)", size=(600,1200))
    push!(P,p)
    vline!([0],c=:black,linestyle=:dash, linewidth=2)
    Plot.save((;xlim="full", key...))
    xlims!(-2,2)
    Plot.save((;xlim=(-2,2), key...))
    xlims!(-1,1)
    Plot.save((;xlim=(-1,1), key...))
end

plot(P..., ylim=(0,1))
shifts = -1:0.025:1

for key in keys(F)
    f = prep(matrixform(F[key]))
    #push_metric!(f, metrics.bitsperspike)
    #push_shiftmetric!(f, best_tau!; metric=:bitsperspike)
    inds = sortperm(f[:bestshift_bitsperspike][:,1])
    bps  = f[:bitsperspike][inds, :]
    XX = hcat([Utils.norm_extrema(b) for b in eachrow(bps)]...)'
    heatmap(shifts, collect(1:size(XX,1)), XX, title="$(key.datacut)", size=(600,1200))
    vline!([0],c=:black,linestyle=:dash, linewidth=2)
    Plot.save((;xlim="full", key...))
    xlims!(-2,2)
    Plot.save((;xlim=(-2,2), key...))
    xlims!(-1,1)
    Plot.save((;xlim=(-1,1), key...))
end
# time: 2022-09-08 15:19:17 EDT
# mode: julia
Plot.setfolder(parent,"fixed_meanmaps2")
# time: 2022-09-08 15:19:17 EDT
# mode: julia
# Visualize the heatmap
# time: 2022-09-08 15:19:17 EDT
# mode: julia
P=[]
# time: 2022-09-08 15:19:17 EDT
# mode: julia
for key in keys(F)
f = prep(matrixform(F[key]))
#push_metric!(f, metrics.bitsperspike)
#push_shiftmetric!(f, best_tau!; metric=:bitsperspike)
inds = sortperm(f[:bestshift_bitsperspike][:,1])
bps  = f[:bitsperspike][inds, :]
XX = hcat([Utils.norm_extrema(b) for b in eachrow(bps)]...)'
p=plot(shifts, vec(mean(XX,dims=1)), title="$(key.datacut)", size=(600,1200))
push!(P,p)
vline!([0],c=:black,linestyle=:dash, linewidth=2)
Plot.save((;xlim="full", key...))
xlims!(-2,2)
Plot.save((;xlim=(-2,2), key...))
    xlims!(-1,1)
        Plot.save((;xlim=(-1,1), key...))
        end
# time: 2022-09-08 15:19:18 EDT
# mode: julia
plot(P..., ylim=(0,1))
# time: 2022-09-08 15:19:18 EDT
# mode: julia
Plot.save("together")
# time: 2022-09-08 15:19:45 EDT
# mode: julia
using Timeshift.shiftmetrics
# time: 2022-09-08 15:19:54 EDT
# mode: julia
Plot.setfolder(parent,"fixed_meanmaps2")
# time: 2022-09-08 15:19:54 EDT
# mode: julia
# Visualize the heatmap
# time: 2022-09-08 15:19:54 EDT
# mode: julia
P=[]
# time: 2022-09-08 15:19:55 EDT
# mode: julia
for key in keys(F)
f = prep(matrixform(F[key]))
#push_metric!(f, metrics.bitsperspike)
#push_shiftmetric!(f, best_tau!; metric=:bitsperspike)
inds = sortperm(f[:bestshift_bitsperspike][:,1])
bps  = f[:bitsperspike][inds, :]
XX = hcat([Utils.norm_extrema(b) for b in eachrow(bps)]...)'
p=plot(shifts, vec(mean(XX,dims=1)), title="$(key.datacut)", size=(600,1200))
push!(P,p)
vline!([0],c=:black,linestyle=:dash, linewidth=2)
Plot.save((;xlim="full", key...))
    xlims!(-2,2)
        Plot.save((;xlim=(-2,2), key...))
            xlims!(-1,1)
                Plot.save((;xlim=(-1,1), key...))
                end
# time: 2022-09-08 15:20:18 EDT
# mode: julia
plot(P..., ylim=(0,1))
# time: 2022-09-08 15:20:18 EDT
# mode: julia
Plot.save("together")
# time: 2022-09-08 15:26:33 EDT
# mode: julia
# Visualize the heatmap
# time: 2022-09-08 15:26:34 EDT
# mode: julia
for key in keys(F)
f = prep(matrixform(F[key]))
#push_metric!(f, metrics.bitsperspike)
#push_shiftmetric!(f, best_tau!; metric=:bitsperspike)
inds = sortperm(f[:bestshift_bitsperspike][:,1])
bps  = f[:bitsperspike][inds, :]
XX = hcat([Utils.norm_extrema(b) for b in eachrow(bps)]...)'
heatmap(shifts, collect(1:size(XX,1)), XX, title="$(key.datacut)", size=(600,1200))
vline!([0],c=:black,linestyle=:dash, linewidth=2)
Plot.save((;xlim="full", key...))
inds = shifts == abs.(shifts) .< 2
heatmap(shifts[inds], collect(1:size(XX,1)), XX[inds,:], title="$(key.datacut)", size=(600,1200))
vline!([0],c=:black,linestyle=:dash, linewidth=2)
    Plot.save((;xlim=(-2,2), key...))
        inds = shifts == abs.(shifts) .< 1
            heatmap(shifts[inds], collect(1:size(XX,1)), XX[inds,:], title="$(key.datacut)", size=(600,1200))
                vline!([0],c=:black,linestyle=:dash, linewidth=2)
                    Plot.save((;xlim=(-1,1), key...))
                    end
# time: 2022-09-08 15:26:46 EDT
# mode: julia
shifts = -2:0.1:2
# time: 2022-09-08 15:26:48 EDT
# mode: julia
# Visualize the heatmap
# time: 2022-09-08 15:26:48 EDT
# mode: julia
for key in keys(F)
f = prep(matrixform(F[key]))
#push_metric!(f, metrics.bitsperspike)
#push_shiftmetric!(f, best_tau!; metric=:bitsperspike)
inds = sortperm(f[:bestshift_bitsperspike][:,1])
bps  = f[:bitsperspike][inds, :]
XX = hcat([Utils.norm_extrema(b) for b in eachrow(bps)]...)'
heatmap(shifts, collect(1:size(XX,1)), XX, title="$(key.datacut)", size=(600,1200))
vline!([0],c=:black,linestyle=:dash, linewidth=2)
Plot.save((;xlim="full", key...))
inds = shifts == abs.(shifts) .< 2
heatmap(shifts[inds], collect(1:size(XX,1)), XX[inds,:], title="$(key.datacut)", size=(600,1200))
vline!([0],c=:black,linestyle=:dash, linewidth=2)
    Plot.save((;xlim=(-2,2), key...))
        inds = shifts == abs.(shifts) .< 1
            heatmap(shifts[inds], collect(1:size(XX,1)), XX[inds,:], title="$(key.datacut)", size=(600,1200))
                vline!([0],c=:black,linestyle=:dash, linewidth=2)
                    Plot.save((;xlim=(-1,1), key...))
                    end
# time: 2022-09-08 15:27:09 EDT
# mode: julia
f = prep(matrixform(F[key]))
# time: 2022-09-08 15:27:17 EDT
# mode: julia
key = first(keys(F))
# time: 2022-09-08 15:27:21 EDT
# mode: julia
f = prep(matrixform(F[key]))
# time: 2022-09-08 15:27:31 EDT
# mode: julia
inds = sortperm(f[:bestshift_bitsperspike][:,1])
# time: 2022-09-08 15:27:32 EDT
# mode: julia
bps  = f[:bitsperspike][inds, :]
# time: 2022-09-08 15:27:33 EDT
# mode: julia
XX = hcat([Utils.norm_extrema(b) for b in eachrow(bps)]...)'
# time: 2022-09-08 15:27:34 EDT
# mode: julia
heatmap(shifts, collect(1:size(XX,1)), XX, title="$(key.datacut)", size=(600,1200))
# time: 2022-09-08 15:27:34 EDT
# mode: julia
vline!([0],c=:black,linestyle=:dash, linewidth=2)
# time: 2022-09-08 15:27:34 EDT
# mode: julia
Plot.save((;xlim="full", key...))
# time: 2022-09-08 15:27:34 EDT
# mode: julia
inds = shifts == abs.(shifts) .< 2
# time: 2022-09-08 15:27:35 EDT
# mode: julia
heatmap(shifts[inds], collect(1:size(XX,1)), XX[inds,:], title="$(key.datacut)", size=(600,1200))
# time: 2022-09-08 15:27:36 EDT
# mode: julia
vline!([0],c=:black,linestyle=:dash, linewidth=2)
# time: 2022-09-08 15:27:36 EDT
# mode: julia
Plot.save((;xlim=(-2,2), key...))
# time: 2022-09-08 15:27:36 EDT
# mode: julia
inds = shifts == abs.(shifts) .< 1
# time: 2022-09-08 15:27:36 EDT
# mode: julia
heatmap(shifts[inds], collect(1:size(XX,1)), XX[inds,:], title="$(key.datacut)", size=(600,1200))
# time: 2022-09-08 15:27:39 EDT
# mode: julia
shifts[inds]
# time: 2022-09-08 15:27:42 EDT
# mode: julia
inds
# time: 2022-09-08 15:27:43 EDT
# mode: julia
inds = shifts == abs.(shifts) .< 1
# time: 2022-09-08 15:27:50 EDT
# mode: julia
shifts
# time: 2022-09-08 15:28:07 EDT
# mode: julia
shift = collect(shifts)
# time: 2022-09-08 15:28:08 EDT
# mode: julia
inds = shifts == abs.(shifts) .< 1
# time: 2022-09-08 15:28:10 EDT
# mode: julia
heatmap(shifts[inds], collect(1:size(XX,1)), XX[inds,:], title="$(key.datacut)", size=(600,1200))
# time: 2022-09-08 15:28:31 EDT
# mode: julia
heatmap(shifts[inds], collect(1:size(XX,1)), XX[:,inds], title="$(key.datacut)", size=(600,1200))
# time: 2022-09-08 15:28:40 EDT
# mode: julia
inds
# time: 2022-09-08 15:28:43 EDT
# mode: julia
size(XX)
# time: 2022-09-08 15:28:46 EDT
# mode: julia
size(inds)
# time: 2022-09-08 15:29:01 EDT
# mode: julia
shifts = collect(shifts)
# time: 2022-09-08 15:30:18 EDT
# mode: julia
f
# time: 2022-09-08 15:30:30 EDT
# mode: julia
f.dims[2]
# time: 2022-09-08 15:30:59 EDT
# mode: julia
shifts = f.dims[2]
# time: 2022-09-08 15:31:04 EDT
# mode: julia
shifts = collect(shifts)
# time: 2022-09-08 15:31:07 EDT
# mode: julia
inds = shifts == abs.(shifts) .< 1
# time: 2022-09-08 15:31:08 EDT
# mode: julia
heatmap(shifts[inds], collect(1:size(XX,1)), XX[inds,:], title="$(key.datacut)", size=(600,1200))
# time: 2022-09-08 15:31:08 EDT
# mode: julia
vline!([0],c=:black,linestyle=:dash, linewidth=2)
# time: 2022-09-08 15:31:10 EDT
# mode: julia
shifts = collect(shifts)
# time: 2022-09-08 15:31:11 EDT
# mode: julia
inds = shifts == abs.(shifts) .< 1
# time: 2022-09-08 15:31:13 EDT
# mode: julia
heatmap(shifts[inds], collect(1:size(XX,1)), XX[inds,:], title="$(key.datacut)", size=(600,1200))
# time: 2022-09-08 15:31:23 EDT
# mode: julia
1:size(XX,1))
# time: 2022-09-08 15:31:34 EDT
# mode: julia
collect(1:size(XX,1))
# time: 2022-09-08 15:31:38 EDT
# mode: julia
shifts[inds]
# time: 2022-09-08 15:31:40 EDT
# mode: julia
inds = shifts == abs.(shifts) .< 1
# time: 2022-09-08 15:31:44 EDT
# mode: julia
any(inds)
# time: 2022-09-08 15:31:47 EDT
# mode: julia
shifsts
# time: 2022-09-08 15:31:49 EDT
# mode: julia
shifts
# time: 2022-09-08 15:31:58 EDT
# mode: julia
inds = abs.(shifts) .< 1
# time: 2022-09-08 15:32:00 EDT
# mode: julia
heatmap(shifts[inds], collect(1:size(XX,1)), XX[inds,:], title="$(key.datacut)", size=(600,1200))
# time: 2022-09-08 15:32:10 EDT
# mode: julia
heatmap(shifts[inds], collect(1:size(XX,1)), XX[:,inds], title="$(key.datacut)", size=(600,1200))
# time: 2022-09-08 15:32:16 EDT
# mode: julia
# Visualize the heatmap
# time: 2022-09-08 15:32:17 EDT
# mode: julia
for key in keys(F)
f = prep(matrixform(F[key]))
#push_metric!(f, metrics.bitsperspike)
#push_shiftmetric!(f, best_tau!; metric=:bitsperspike)
shifts = f.dims[2]
inds = sortperm(f[:bestshift_bitsperspike][:,1])
bps  = f[:bitsperspike][inds, :]
XX = hcat([Utils.norm_extrema(b) for b in eachrow(bps)]...)'
heatmap(shifts, collect(1:size(XX,1)), XX, title="$(key.datacut)", size=(600,1200))
vline!([0],c=:black,linestyle=:dash, linewidth=2)
Plot.save((;xlim="full", key...))
inds = shifts == abs.(shifts) .< 2
heatmap(shifts[inds], collect(1:size(XX,1)), XX[inds,:], title="$(key.datacut)", size=(600,1200))
vline!([0],c=:black,linestyle=:dash, linewidth=2)
Plot.save((;xlim=(-2,2), key...))
    shifts = collect(shifts)
        inds = abs.(shifts) .< 1
            heatmap(shifts[inds], collect(1:size(XX,1)), XX[:,inds], title="$(key.datacut)", size=(600,1200))
                vline!([0],c=:black,linestyle=:dash, linewidth=2)
                    Plot.save((;xlim=(-1,1), key...))
                    end
# time: 2022-09-08 15:32:53 EDT
# mode: julia
# Visualize the heatmap
# time: 2022-09-08 15:32:54 EDT
# mode: julia
for key in keys(F)
f = prep(matrixform(F[key]))
#push_metric!(f, metrics.bitsperspike)
#push_shiftmetric!(f, best_tau!; metric=:bitsperspike)
shifts = vec(f.dims[2])
inds = sortperm(f[:bestshift_bitsperspike][:,1])
bps  = f[:bitsperspike][inds, :]
XX = hcat([Utils.norm_extrema(b) for b in eachrow(bps)]...)'
heatmap(shifts, collect(1:size(XX,1)), XX, title="$(key.datacut)", size=(600,1200))
vline!([0],c=:black,linestyle=:dash, linewidth=2)
Plot.save((;xlim="full", key...))
    inds = shifts == abs.(shifts) .< 2
        heatmap(shifts[inds], collect(1:size(XX,1)), XX[inds,:], title="$(key.datacut)", size=(600,1200))
            vline!([0],c=:black,linestyle=:dash, linewidth=2)
                Plot.save((;xlim=(-2,2), key...))
                    shifts = collect(shifts)
                        inds = abs.(shifts) .< 1
                            heatmap(shifts[inds], collect(1:size(XX,1)), XX[:,inds], title="$(key.datacut)", size=(600,1200))
                                vline!([0],c=:black,linestyle=:dash, linewidth=2)
                                    Plot.save((;xlim=(-1,1), key...))
                                    end
# time: 2022-09-08 15:33:15 EDT
# mode: julia
f.dims[2]
# time: 2022-09-08 15:33:20 EDT
# mode: julia
f.dims[2].val
# time: 2022-09-08 15:33:25 EDT
# mode: julia
vec(f.dims[2])
# time: 2022-09-08 15:33:30 EDT
# mode: julia
f.dims[2]
# time: 2022-09-08 15:33:33 EDT
# mode: julia
d=f.dims[2]
# time: 2022-09-08 15:33:39 EDT
# mode: julia
d=f.dims[2].val
# time: 2022-09-08 15:33:42 EDT
# mode: julia
d=vec(f.dims[2].val)
# time: 2022-09-08 15:34:10 EDT
# mode: julia
f.dims[2].val
# time: 2022-09-08 15:34:12 EDT
# mode: julia
f.dims[2].val.val
# time: 2022-09-08 15:34:20 EDT
# mode: julia
Vector(f.dims[2])
# time: 2022-09-08 15:34:26 EDT
# mode: julia
f.dims[2]
# time: 2022-09-08 15:34:30 EDT
# mode: julia
f.dims[2].val
# time: 2022-09-08 15:34:35 EDT
# mode: julia
v=f.dims[2].val
# time: 2022-09-08 15:34:38 EDT
# mode: julia
v=f.dims[2].val.data
# time: 2022-09-08 15:34:50 EDT
# mode: julia
# Visualize the heatmap
# time: 2022-09-08 15:34:51 EDT
# mode: julia
for key in keys(F)
f = prep(matrixform(F[key]))
#push_metric!(f, metrics.bitsperspike)
#push_shiftmetric!(f, best_tau!; metric=:bitsperspike)
shifts = vec(f.dims[2].val.data)
inds = sortperm(f[:bestshift_bitsperspike][:,1])
bps  = f[:bitsperspike][inds, :]
XX = hcat([Utils.norm_extrema(b) for b in eachrow(bps)]...)'
heatmap(shifts, collect(1:size(XX,1)), XX, title="$(key.datacut)", size=(600,1200))
vline!([0],c=:black,linestyle=:dash, linewidth=2)
Plot.save((;xlim="full", key...))
inds = shifts == abs.(shifts) .< 2
heatmap(shifts[inds], collect(1:size(XX,1)), XX[inds,:], title="$(key.datacut)", size=(600,1200))
    vline!([0],c=:black,linestyle=:dash, linewidth=2)
        Plot.save((;xlim=(-2,2), key...))
            shifts = collect(shifts)
                inds = abs.(shifts) .< 1
                    heatmap(shifts[inds], collect(1:size(XX,1)), XX[:,inds], title="$(key.datacut)", size=(600,1200))
                        vline!([0],c=:black,linestyle=:dash, linewidth=2)
                            Plot.save((;xlim=(-1,1), key...))
                            end
# time: 2022-09-08 15:35:23 EDT
# mode: julia
shifts = collect(shifts)
# time: 2022-09-08 15:35:25 EDT
# mode: julia
inds = abs.(shifts) .< 1
# time: 2022-09-08 15:35:26 EDT
# mode: julia
heatmap(shifts[inds], collect(1:size(XX,1)), XX[:,inds], title="$(key.datacut)", size=(600,1200))
# time: 2022-09-08 15:35:42 EDT
# mode: julia
inds = abs.(shifts) .< 2
# time: 2022-09-08 15:35:42 EDT
# mode: julia
heatmap(shifts[inds], collect(1:size(XX,1)), XX[inds,:], title="$(key.datacut)", size=(600,1200))
# time: 2022-09-08 15:35:42 EDT
# mode: julia
vline!([0],c=:black,linestyle=:dash, linewidth=2)
# time: 2022-09-08 15:35:42 EDT
# mode: julia
Plot.save((;xlim=(-2,2), key...))
# time: 2022-09-08 15:35:58 EDT
# mode: julia
inds = abs.(shifts) .< 2
# time: 2022-09-08 15:35:58 EDT
# mode: julia
heatmap(shifts[inds], collect(1:size(XX,1)), XX[:,inds], title="$(key.datacut)", size=(600,1200))
# time: 2022-09-08 15:35:58 EDT
# mode: julia
vline!([0],c=:black,linestyle=:dash, linewidth=2)
# time: 2022-09-08 15:36:03 EDT
# mode: julia
# Visualize the heatmap
# time: 2022-09-08 15:36:04 EDT
# mode: julia
for key in keys(F)
f = prep(matrixform(F[key]))
#push_metric!(f, metrics.bitsperspike)
#push_shiftmetric!(f, best_tau!; metric=:bitsperspike)
shifts = vec(f.dims[2].val.data)
inds = sortperm(f[:bestshift_bitsperspike][:,1])
bps  = f[:bitsperspike][inds, :]
XX = hcat([Utils.norm_extrema(b) for b in eachrow(bps)]...)'
heatmap(shifts, collect(1:size(XX,1)), XX, title="$(key.datacut)", size=(600,1200))
vline!([0],c=:black,linestyle=:dash, linewidth=2)
Plot.save((;xlim="full", key...))

    inds = abs.(shifts) .< 2
        heatmap(shifts[inds], collect(1:size(XX,1)), XX[:,inds], title="$(key.datacut)", size=(600,1200))
            vline!([0],c=:black,linestyle=:dash, linewidth=2)
                Plot.save((;xlim=(-2,2), key...))
                
                    inds = abs.(shifts) .< 1
                        heatmap(shifts[inds], collect(1:size(XX,1)), XX[:,inds], title="$(key.datacut)", size=(600,1200))
                            vline!([0],c=:black,linestyle=:dash, linewidth=2)
                                Plot.save((;xlim=(-1,1), key...))
                                end
# time: 2022-09-08 15:49:30 EDT
# mode: julia
@info "keys" k1 k2
# time: 2022-09-08 15:49:32 EDT
# mode: julia
key1,key2 = keyz[k1],keyz[k2]
# time: 2022-09-08 15:49:32 EDT
# mode: julia
f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]))
# time: 2022-09-08 15:49:32 EDT
# mode: julia
u = intersect(f1.dims[1], f2.dims[1])
# time: 2022-09-08 15:49:36 EDT
# mode: julia
key1,key2 = keyz[k1],keyz[k2]
# time: 2022-09-08 15:49:36 EDT
# mode: julia
f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]))
# time: 2022-09-08 15:49:36 EDT
# mode: julia
u = intersect(f1.dims[1], f2.dims[1])
# time: 2022-09-08 15:49:38 EDT
# mode: julia
@info "keys" k1 k2
# time: 2022-09-08 15:49:46 EDT
# mode: julia
k1, k2 = :cue, :memory
# time: 2022-09-08 15:49:52 EDT
# mode: julia
key1,key2 = keyz[k1],keyz[k2]
# time: 2022-09-08 15:49:55 EDT
# mode: julia
@info "keys" k1 k2
# time: 2022-09-08 15:49:56 EDT
# mode: julia
key1,key2 = keyz[k1],keyz[k2]
# time: 2022-09-08 15:50:00 EDT
# mode: julia
keyz
# time: 2022-09-08 15:50:04 EDT
# mode: julia
keyz = Dict(key.datacut=>key for key in keys(F))
# time: 2022-09-08 15:50:08 EDT
# mode: julia
key1,key2 = keyz[k1],keyz[k2]
# time: 2022-09-08 15:50:09 EDT
# mode: julia
f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]))
# time: 2022-09-08 15:50:18 EDT
# mode: julia
u = intersect(f1.dims[1], f2.dims[1])
# time: 2022-09-08 15:50:20 EDT
# mode: julia
vals = f1[unit=At(u)][:bestshift_bitsperspike],
f2[unit=At(u)][:bestshift_bitsperspike]
# time: 2022-09-08 15:50:30 EDT
# mode: julia
using DimensionalData
# time: 2022-09-08 15:50:35 EDT
# mode: julia
vals = f1[unit=At(u)][:bestshift_bitsperspike],
f2[unit=At(u)][:bestshift_bitsperspike]
# time: 2022-09-08 15:50:35 EDT
# mode: julia
# HISTOGRAM ("PDF")
# time: 2022-09-08 15:50:38 EDT
# mode: julia
plot_histdist = 
histogram(shiftdist; bins=edges,  xlim, normalize=:probability, alpha=0.5, label="", xlabel="shift", ylabel="pdf")
# time: 2022-09-08 15:50:38 EDT
# mode: julia
vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
# time: 2022-09-08 15:51:02 EDT
# mode: julia
diff(vals)
# time: 2022-09-08 15:51:16 EDT
# mode: julia
shifdist = vals[2] - vals[1]
# time: 2022-09-08 15:51:33 EDT
# mode: julia
shifdist = (vals[2] - vals[1])[:,1]
# time: 2022-09-08 15:51:41 EDT
# mode: julia
plot_histdist = 
histogram(shiftdist; bins=edges,  xlim, normalize=:probability, alpha=0.5, label="", xlabel="shift", ylabel="pdf")
# time: 2022-09-08 15:52:25 EDT
# mode: julia
edges = LinRange(minimum(f1[:shift]), maximum(f1[:shift]), bins+1)
# time: 2022-09-08 15:52:41 EDT
# mode: julia
bins = 25
# time: 2022-09-08 15:52:45 EDT
# mode: julia
edges = LinRange(minimum(f1[:shift]), maximum(f1[:shift]), bins+1)
# time: 2022-09-08 15:52:48 EDT
# mode: julia
shifdist = (vals[2] - vals[1])[:,1]
# time: 2022-09-08 15:52:56 EDT
# mode: julia
plot_histdist = 
histogram(shiftdist; bins=edges,  xlim, normalize=:probability, alpha=0.5, label="", xlabel="shift", ylabel="pdf")
# time: 2022-09-08 15:52:56 EDT
# mode: julia
vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
# time: 2022-09-08 15:53:20 EDT
# mode: julia
xlim = (edges[begin]. edges[end])
# time: 2022-09-08 15:53:23 EDT
# mode: julia
edges = LinRange(minimum(f1[:shift]), maximum(f1[:shift]), bins+1)
# time: 2022-09-08 15:53:25 EDT
# mode: julia
xlim = (edges[begin]. edges[end])
# time: 2022-09-08 15:53:32 EDT
# mode: julia
xlim = (edges[begin], edges[end])
# time: 2022-09-08 15:53:35 EDT
# mode: julia
shifdist = (vals[2] - vals[1])[:,1]
# time: 2022-09-08 15:53:35 EDT
# mode: julia
# HISTOGRAM ("PDF")
# time: 2022-09-08 15:53:35 EDT
# mode: julia
plot_histdist = 
histogram(shiftdist; bins=edges,  xlim, normalize=:probability, alpha=0.5, label="", xlabel="shift", ylabel="pdf")
# time: 2022-09-08 15:53:35 EDT
# mode: julia
vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
# time: 2022-09-08 15:53:38 EDT
# mode: julia
shifdist = (vals[2] - vals[1])[:,1]
# time: 2022-09-08 15:53:43 EDT
# mode: julia
shiftdist = (vals[2] - vals[1])[:,1]
# time: 2022-09-08 15:53:45 EDT
# mode: julia
plot_histdist = 
histogram(shiftdist; bins=edges,  xlim, normalize=:probability, alpha=0.5, label="", xlabel="shift", ylabel="pdf")
# time: 2022-09-08 15:53:49 EDT
# mode: julia
vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
# time: 2022-09-08 15:53:54 EDT
# mode: julia
dist = fit(Histogram, shiftdist, edges)
# time: 2022-09-08 15:53:54 EDT
# mode: julia
dist = StatsBase.normalize(dist, mode=:probability)
# time: 2022-09-08 15:53:54 EDT
# mode: julia
plot_histdist_cdf = 
bar(dist.edges, cumsum(dist.weights); xlim, alpha=0.5, label="", xlabel="shift", ylabel="cdf")
# time: 2022-09-08 15:53:54 EDT
# mode: julia
vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
# time: 2022-09-08 15:54:05 EDT
# mode: julia
# Visualize the heatmap
# time: 2022-09-08 15:54:05 EDT
# mode: julia
keyz = Dict(key.datacut=>key for key in keys(F))
# time: 2022-09-08 15:54:05 EDT
# mode: julia
Plot.setfolder(parent,"hist_diff_correct")
# time: 2022-09-08 15:54:05 EDT
# mode: julia
bins = 25
# time: 2022-09-08 15:54:07 EDT
# mode: julia
for (k1, k2) in zip((:nontask, :cue), (:task,:memory))
@info "keys" k1 k2
key1,key2 = keyz[k1],keyz[k2]
f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]))
u = intersect(f1.dims[1], f2.dims[1])

    vals = f1[unit=At(u)][:bestshift_bitsperspike],
               f2[unit=At(u)][:bestshift_bitsperspike]
               
                  edges = LinRange(minimum(f1[:shift]), maximum(f1[:shift]), bins+1)
                     xlim = (edges[begin], edges[end])
                     
                         shiftdist = (vals[2] - vals[1])[:,1]
                             # HISTOGRAM ("PDF")
                                 plot_histdist = 
                                         histogram(shiftdist; bins=edges,  xlim, normalize=:probability, alpha=0.5, label="", xlabel="shift", ylabel="pdf")
                                             vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
                                             
                                                 # HISTOGRAM ("CDF")
                                                     dist = fit(Histogram, shiftdist, edges)
                                                         dist = StatsBase.normalize(dist, mode=:probability)
                                                             plot_histdist_cdf = 
                                                                     bar(dist.edges, cumsum(dist.weights); xlim, alpha=0.5, label="", xlabel="shift", ylabel="cdf")
                                                                         vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
                                                                             #scatter(vec.(vals)...)
                                                                                 #plot!(-2:2,-2:2,c=:white)
                                                                                     Plot.save((;k1,k2))
                                                                                     end
# time: 2022-09-08 15:54:24 EDT
# mode: julia
# Visualize the heatmap
# time: 2022-09-08 15:54:24 EDT
# mode: julia
using StatsBase
# time: 2022-09-08 15:54:24 EDT
# mode: julia
keyz = Dict(key.datacut=>key for key in keys(F))
# time: 2022-09-08 15:54:24 EDT
# mode: julia
Plot.setfolder(parent,"hist_diff_correct")
# time: 2022-09-08 15:54:24 EDT
# mode: julia
bins = 25
# time: 2022-09-08 15:54:25 EDT
# mode: julia
for (k1, k2) in zip((:nontask, :cue), (:task,:memory))
@info "keys" k1 k2
key1,key2 = keyz[k1],keyz[k2]
f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]))
u = intersect(f1.dims[1], f2.dims[1])

    vals = f1[unit=At(u)][:bestshift_bitsperspike],
               f2[unit=At(u)][:bestshift_bitsperspike]
               
                  edges = LinRange(minimum(f1[:shift]), maximum(f1[:shift]), bins+1)
                     xlim = (edges[begin], edges[end])
                     
                         shiftdist = (vals[2] - vals[1])[:,1]
                             # HISTOGRAM ("PDF")
                                 plot_histdist = 
                                         histogram(shiftdist; bins=edges,  xlim, normalize=:probability, alpha=0.5, label="", xlabel="shift", ylabel="pdf")
                                             vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
                                             
                                                 # HISTOGRAM ("CDF")
                                                     dist = fit(Histogram, shiftdist, edges)
                                                         dist = StatsBase.normalize(dist, mode=:probability)
                                                             plot_histdist_cdf = 
                                                                     bar(dist.edges, cumsum(dist.weights); xlim, alpha=0.5, label="", xlabel="shift", ylabel="cdf")
                                                                         vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
                                                                             #scatter(vec.(vals)...)
                                                                                 #plot!(-2:2,-2:2,c=:white)
                                                                                     Plot.save((;k1,k2))
                                                                                     end
# time: 2022-09-08 15:55:15 EDT
# mode: julia
# Visualize the heatmap
# time: 2022-09-08 15:55:15 EDT
# mode: julia
using StatsBase
# time: 2022-09-08 15:55:15 EDT
# mode: julia
keyz = Dict(key.datacut=>key for key in keys(F))
# time: 2022-09-08 15:55:15 EDT
# mode: julia
Plot.setfolder(parent,"hist_diff_correct")
# time: 2022-09-08 15:55:15 EDT
# mode: julia
bins = 25
# time: 2022-09-08 15:55:17 EDT
# mode: julia
for (k1, k2) in zip((:nontask, :cue), (:task,:memory))
@info "keys" k1 k2
key1,key2 = keyz[k1],keyz[k2]
f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]))
u = intersect(f1.dims[1], f2.dims[1])

    vals = f1[unit=At(u)][:bestshift_bitsperspike],
               f2[unit=At(u)][:bestshift_bitsperspike]
               
                  edges = LinRange(minimum(f1[:shift]), maximum(f1[:shift]), bins+1)
                     xlim = (edges[begin], edges[end])
                     
                         shiftdist = (vals[2] - vals[1])[:,1]
                             # HISTOGRAM ("PDF")
                                 plot_histdist = 
                                         histogram(shiftdist; bins=edges,  xlim, normalize=:probability, alpha=0.5, label="", xlabel="shift", ylabel="pdf")
                                             vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
                                             
                                                 # HISTOGRAM ("CDF")
                                                     dist = fit(Histogram, shiftdist, edges)
                                                         dist = StatsBase.normalize(dist, mode=:probability)
                                                             plot_histdist_cdf = 
                                                                     bar(dist.edges, cumsum(dist.weights); xlim, alpha=0.5, label="", xlabel="shift", ylabel="cdf")
                                                                         vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
                                                                             plot(plot_histdist, plot_histdist_cdf, layout=grid(2,1))
                                                                                 #scatter(vec.(vals)...)
                                                                                     #plot!(-2:2,-2:2,c=:white)
                                                                                         Plot.save((;k1,k2))
                                                                                         end
# time: 2022-09-08 16:01:09 EDT
# mode: julia
future = shifdist .> 0
# time: 2022-09-08 16:01:09 EDT
# mode: julia
past = shifdist .< 0
# time: 2022-09-08 16:01:26 EDT
# mode: julia
future = sum(shifdist .> 0)
# time: 2022-09-08 16:01:26 EDT
# mode: julia
past = sum(shifdist .< 0)
# time: 2022-09-08 16:01:58 EDT
# mode: julia
shiftdist = (vals[2] - vals[1])[:,1]
# time: 2022-09-08 16:01:58 EDT
# mode: julia
# HISTOGRAM ("PDF")
# time: 2022-09-08 16:01:58 EDT
# mode: julia
plot_histdist = 
histogram(shiftdist; bins=edges,  xlim, normalize=:probability, alpha=0.5, label="", xlabel="shift", ylabel="pdf")
# time: 2022-09-08 16:01:58 EDT
# mode: julia
vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
# time: 2022-09-08 16:01:58 EDT
# mode: julia
# HISTOGRAM ("CDF")
# time: 2022-09-08 16:01:58 EDT
# mode: julia
dist = fit(Histogram, shiftdist, edges)
# time: 2022-09-08 16:01:58 EDT
# mode: julia
dist = StatsBase.normalize(dist, mode=:probability)
# time: 2022-09-08 16:01:58 EDT
# mode: julia
plot_histdist_cdf = 
bar(dist.edges, cumsum(dist.weights); xlim, alpha=0.5, label="", xlabel="shift", ylabel="cdf")
# time: 2022-09-08 16:01:58 EDT
# mode: julia
vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
# time: 2022-09-08 16:01:58 EDT
# mode: julia
future = sum(shifdist .> 0)
# time: 2022-09-08 16:01:58 EDT
# mode: julia
past = sum(shifdist .< 0)
# time: 2022-09-08 16:01:58 EDT
# mode: julia
plot(plot_histdist, plot_histdist_cdf, layout=grid(2,1))
# time: 2022-09-08 16:03:16 EDT
# mode: julia
k1=:nontask, k2=:task
# time: 2022-09-08 16:03:22 EDT
# mode: julia
k1=:nontask; k2=:task
# time: 2022-09-08 16:03:54 EDT
# mode: julia
key1,key2 = keyz[k1],keyz[k2]
# time: 2022-09-08 16:03:54 EDT
# mode: julia
f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]))
# time: 2022-09-08 16:03:58 EDT
# mode: julia
u = intersect(f1.dims[1], f2.dims[1])
# time: 2022-09-08 16:03:58 EDT
# mode: julia
vals = f1[unit=At(u)][:bestshift_bitsperspike],
f2[unit=At(u)][:bestshift_bitsperspike]
# time: 2022-09-08 16:03:58 EDT
# mode: julia
edges = LinRange(minimum(f1[:shift]), maximum(f1[:shift]), bins+1)
# time: 2022-09-08 16:03:58 EDT
# mode: julia
xlim = (edges[begin], edges[end])
# time: 2022-09-08 16:03:58 EDT
# mode: julia
shiftdist = (vals[2] - vals[1])[:,1]
# time: 2022-09-08 16:03:58 EDT
# mode: julia
# HISTOGRAM ("PDF")
# time: 2022-09-08 16:03:58 EDT
# mode: julia
plot_histdist = 
histogram(shiftdist; bins=edges,  xlim, normalize=:probability, alpha=0.5, label="", xlabel="shift", ylabel="pdf")
# time: 2022-09-08 16:03:58 EDT
# mode: julia
vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
# time: 2022-09-08 16:03:58 EDT
# mode: julia
# HISTOGRAM ("CDF")
# time: 2022-09-08 16:03:58 EDT
# mode: julia
dist = fit(Histogram, shiftdist, edges)
# time: 2022-09-08 16:03:58 EDT
# mode: julia
dist = StatsBase.normalize(dist, mode=:probability)
# time: 2022-09-08 16:03:58 EDT
# mode: julia
plot_histdist_cdf = 
bar(dist.edges, cumsum(dist.weights); xlim, alpha=0.5, label="", xlabel="shift", ylabel="cdf")
# time: 2022-09-08 16:03:58 EDT
# mode: julia
vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
# time: 2022-09-08 16:03:58 EDT
# mode: julia
future = sum(shifdist .> 0)
# time: 2022-09-08 16:03:58 EDT
# mode: julia
past = sum(shifdist .< 0)
# time: 2022-09-08 16:04:44 EDT
# mode: julia
f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]));
# time: 2022-09-08 16:04:47 EDT
# mode: julia
u = intersect(f1.dims[1], f2.dims[1])
# time: 2022-09-08 16:04:48 EDT
# mode: julia
u
# time: 2022-09-08 16:04:52 EDT
# mode: julia
k
# time: 2022-09-08 16:04:53 EDT
# mode: julia
k1
# time: 2022-09-08 16:04:54 EDT
# mode: julia
k2
# time: 2022-09-08 16:04:59 EDT
# mode: julia
edges = LinRange(minimum(f1[:shift]), maximum(f1[:shift]), bins+1)
# time: 2022-09-08 16:04:59 EDT
# mode: julia
xlim = (edges[begin], edges[end])
# time: 2022-09-08 16:05:06 EDT
# mode: julia
shiftdist = (vals[2] - vals[1])[:,1]
# time: 2022-09-08 16:05:06 EDT
# mode: julia
# HISTOGRAM ("PDF")
# time: 2022-09-08 16:05:06 EDT
# mode: julia
plot_histdist = 
histogram(shiftdist; bins=edges,  xlim, normalize=:probability, alpha=0.5, label="", xlabel="shift", ylabel="pdf")
# time: 2022-09-08 16:05:06 EDT
# mode: julia
vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
# time: 2022-09-08 16:05:06 EDT
# mode: julia
# HISTOGRAM ("CDF")
# time: 2022-09-08 16:05:06 EDT
# mode: julia
dist = fit(Histogram, shiftdist, edges)
# time: 2022-09-08 16:05:06 EDT
# mode: julia
dist = StatsBase.normalize(dist, mode=:probability)
# time: 2022-09-08 16:05:06 EDT
# mode: julia
plot_histdist_cdf = 
bar(dist.edges, cumsum(dist.weights); xlim, alpha=0.5, label="", xlabel="shift", ylabel="cdf")
# time: 2022-09-08 16:05:06 EDT
# mode: julia
vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
# time: 2022-09-08 16:05:06 EDT
# mode: julia
future = sum(shifdist .> 0)
# time: 2022-09-08 16:05:06 EDT
# mode: julia
past = sum(shifdist .< 0)
# time: 2022-09-08 16:05:06 EDT
# mode: julia
plot(plot_histdist, plot_histdist_cdf, layout=grid(2,1))
# time: 2022-09-08 16:05:14 EDT
# mode: julia
k1,k2=:cue,:mem
# time: 2022-09-08 16:05:16 EDT
# mode: julia
k1,k2=:cue,:memory
# time: 2022-09-08 16:05:22 EDT
# mode: julia
@info "keys" k1 k2
# time: 2022-09-08 16:05:22 EDT
# mode: julia
key1,key2 = keyz[k1],keyz[k2]
# time: 2022-09-08 16:05:22 EDT
# mode: julia
f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]));
# time: 2022-09-08 16:05:25 EDT
# mode: julia
u = intersect(f1.dims[1], f2.dims[1])
# time: 2022-09-08 16:05:25 EDT
# mode: julia
vals = f1[unit=At(u)][:bestshift_bitsperspike],
f2[unit=At(u)][:bestshift_bitsperspike]
# time: 2022-09-08 16:05:25 EDT
# mode: julia
edges = LinRange(minimum(f1[:shift]), maximum(f1[:shift]), bins+1)
# time: 2022-09-08 16:05:25 EDT
# mode: julia
xlim = (edges[begin], edges[end])
# time: 2022-09-08 16:05:25 EDT
# mode: julia
shiftdist = (vals[2] - vals[1])[:,1]
# time: 2022-09-08 16:05:25 EDT
# mode: julia
# HISTOGRAM ("PDF")
# time: 2022-09-08 16:05:25 EDT
# mode: julia
plot_histdist = 
histogram(shiftdist; bins=edges,  xlim, normalize=:probability, alpha=0.5, label="", xlabel="shift", ylabel="pdf")
# time: 2022-09-08 16:05:25 EDT
# mode: julia
vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
# time: 2022-09-08 16:05:25 EDT
# mode: julia
# HISTOGRAM ("CDF")
# time: 2022-09-08 16:05:25 EDT
# mode: julia
dist = fit(Histogram, shiftdist, edges)
# time: 2022-09-08 16:05:25 EDT
# mode: julia
dist = StatsBase.normalize(dist, mode=:probability)
# time: 2022-09-08 16:05:25 EDT
# mode: julia
plot_histdist_cdf = 
bar(dist.edges, cumsum(dist.weights); xlim, alpha=0.5, label="", xlabel="shift", ylabel="cdf")
# time: 2022-09-08 16:05:25 EDT
# mode: julia
vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
# time: 2022-09-08 16:05:25 EDT
# mode: julia
future = sum(shifdist .> 0)
# time: 2022-09-08 16:05:25 EDT
# mode: julia
past = sum(shifdist .< 0)
# time: 2022-09-08 16:05:41 EDT
# mode: julia
future = sum(shiftdist .> 0)
# time: 2022-09-08 16:05:41 EDT
# mode: julia
past = sum(shiftdist .< 0)
# time: 2022-09-08 16:05:47 EDT
# mode: julia
k1=:nontask; k2=:task
# time: 2022-09-08 16:05:54 EDT
# mode: julia
@info "keys" k1 k2
# time: 2022-09-08 16:05:54 EDT
# mode: julia
key1,key2 = keyz[k1],keyz[k2]
# time: 2022-09-08 16:05:54 EDT
# mode: julia
f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]));
# time: 2022-09-08 16:05:57 EDT
# mode: julia
u = intersect(f1.dims[1], f2.dims[1])
# time: 2022-09-08 16:05:57 EDT
# mode: julia
vals = f1[unit=At(u)][:bestshift_bitsperspike],
f2[unit=At(u)][:bestshift_bitsperspike]
# time: 2022-09-08 16:05:57 EDT
# mode: julia
edges = LinRange(minimum(f1[:shift]), maximum(f1[:shift]), bins+1)
# time: 2022-09-08 16:05:57 EDT
# mode: julia
xlim = (edges[begin], edges[end])
# time: 2022-09-08 16:05:57 EDT
# mode: julia
shiftdist = (vals[2] - vals[1])[:,1]
# time: 2022-09-08 16:05:57 EDT
# mode: julia
# HISTOGRAM ("PDF")
# time: 2022-09-08 16:05:57 EDT
# mode: julia
plot_histdist = 
histogram(shiftdist; bins=edges,  xlim, normalize=:probability, alpha=0.5, label="", xlabel="shift", ylabel="pdf")
# time: 2022-09-08 16:05:57 EDT
# mode: julia
vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
# time: 2022-09-08 16:05:57 EDT
# mode: julia
# HISTOGRAM ("CDF")
# time: 2022-09-08 16:05:57 EDT
# mode: julia
dist = fit(Histogram, shiftdist, edges)
# time: 2022-09-08 16:05:57 EDT
# mode: julia
dist = StatsBase.normalize(dist, mode=:probability)
# time: 2022-09-08 16:05:57 EDT
# mode: julia
plot_histdist_cdf = 
bar(dist.edges, cumsum(dist.weights); xlim, alpha=0.5, label="", xlabel="shift", ylabel="cdf")
# time: 2022-09-08 16:05:57 EDT
# mode: julia
vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
# time: 2022-09-08 16:05:57 EDT
# mode: julia
future = sum(shiftdist .> 0)
# time: 2022-09-08 16:05:57 EDT
# mode: julia
past = sum(shiftdist .< 0)
# time: 2022-09-08 16:05:57 EDT
# mode: julia
plot(plot_histdist, plot_histdist_cdf, layout=grid(2,1))
# time: 2022-09-08 16:12:05 EDT
# mode: julia
using SoftGlobalScope
# time: 2022-09-08 16:12:26 EDT
# mode: julia
@softscope for  w in [8]
Plot.setfolder("timeshift","xyG-xywidth=$w")

    widths["stopWell"] = 0.50
widths["x"] = w
widths["y"] = w
Utils.filtreg.register(beh,spikes,on="time",transfer=["stopWell"])

    # NOTE :
        # IF you want to use :nonatsk in here, you need to have pseudo-stopWell
            # determined. This means using state-transition history to determine
                # the best fiting route @ a given time
# time: 2022-09-08 16:12:27 EDT
# mode: julia
#@showprogress for datacut in [:all, :task, :cue, :memory]
# time: 2022-09-08 16:12:36 EDT
# mode: julia
@showprogress for datacut in [:all, :task, :cue, :memory, :correct, :error]
        shifted = Timeshift.shifted_fields(dropmissing(beh,:stopWell),
                                               dropmissing(spikes, :stopWell), 
                                                                                      shifts, ["x","y"];
                                                                                                                             shiftbeh=false,
                                                                                                                                                                    widths, 
                                                                                                                                                                                                           adaptive=false,
                                                                                                                                                                                                                                                  metricfuncs=[metrics.bitsperspike,metrics.totalcount,metrics.maxrate,metrics.meanrate],
                                                                                                                                                                                                                                                                                         filters=filts[datacut], 
                                                                                                                                                                                                                                                                                                                                thresh);
                                                                                                                                                                                                                                                                                                                                        shifted_wg = Timeshift.shifted_fields(dropmissing(beh,:stopWell),
                                                                                                                                                                                                                                                                                                                                                                               dropmissing(spikes, :stopWell), 
                                                                                                                                                                                                                                                                                                                                                                                                                      shifts, ["x","y","stopWell"];
                                                                                                                                                                                                                                                                                                                                                                                                                                                             shiftbeh=false,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    widths, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           adaptive=false,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  metricfuncs=[metrics.bitsperspike,metrics.totalcount,metrics.maxrate,metrics.meanrate],
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         filters=filts[datacut], 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                thresh);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        fwg, f = matrixform(ShiftedFields(shifted_wg)), matrixform(ShiftedFields(shifted))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                serialize(datadir("exp_pro", "xyG-$datacut-$w-fmat"), (;f, fwg))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        sh = collect(shifts)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                #f, fwg = deserialize(datadir("exp_pro", "xyG-$datacut-$w-fmat"))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        push_shiftmetric!(fwg, best_tau!; metric=:bitsperspike)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                gd = [Utils.squeeze(nanmean(nanmean(f.rate; dims=1), dims=2))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              for f in fwg]
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      gm = [nanmaximum(g)-nanminimum(g) for g in gd]
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              g_count = [Utils.squeeze(nansum(nansum(f.count; dims=1), dims=2))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  for f in fwg]
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          g_occ_count = [Utils.squeeze(nansum(nansum(f.occ.count; dims=1), dims=2))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              for f in fwg]
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      gc = [nanmaximum(c./o)-nanminimum(c./o)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    for (c,o) in zip(g_count, g_occ_count)]
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            fwg[:bestshift_bitsperspike][fwg[:bestshift_bitsperspike] .< 0.5] .= NaN
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    gm = DimArray(gm, fwg.dims)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            B = []
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    for (r_fwg, r_gm) in zip(eachrow(fwg),eachrow(gm))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                samp = [r_fwg[shift=At(r_fwg[:bestshift_bitsperspike][1])][:bestshift_bitsperspike],
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             r_gm[shift=At(r_fwg[:bestshift_bitsperspike][1])]]
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         push!(B, samp)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 end
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         B = hcat(B...)'
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 B[:,1] = abs.(B[:,1])
# time: 2022-09-08 16:12:36 EDT
# mode: julia
    se(x)  = std(x)/sqrt(length(x))
# time: 2022-09-08 16:12:36 EDT
# mode: julia
    inds = B[:,1] .> 0.5
# time: 2022-09-08 16:12:36 EDT
# mode: julia
    chunks = [(B[(!).(inds),2]), (B[inds,2])]
# time: 2022-09-08 16:12:36 EDT
# mode: julia
    t = OneWayANOVATest(chunks...)
# time: 2022-09-08 16:12:36 EDT
# mode: julia
    bar(1:2, mean.(chunks), yerror=se.(chunks), title="pval=$(pvalue(t))")
# time: 2022-09-08 16:12:36 EDT
# mode: julia
    Plot.save("summary, datacut=$datacut, w=$w, pval=$(pvalue(t))")
# time: 2022-09-08 16:12:36 EDT
# mode: julia
    inds = B[:,1] .> 1
# time: 2022-09-08 16:12:36 EDT
# mode: julia
    chunks = [(B[(!).(inds),2]), (B[inds,2])]
# time: 2022-09-08 16:12:36 EDT
# mode: julia
    t = OneWayANOVATest(chunks...)
# time: 2022-09-08 16:12:36 EDT
# mode: julia
    bar(1:2, mean.(chunks), yerror=se.(chunks), title="pval=$(pvalue(t))")
# time: 2022-09-08 16:12:36 EDT
# mode: julia
    Plot.save("summary, <1, datacut=$datacut, w=$w, pval=$(pvalue(t))")
# time: 2022-09-08 16:12:37 EDT
# mode: julia
    ##scatter(vec(fwg[:bestshift_bitsperspike]), vec(gm), title="$datacut")
# time: 2022-09-08 16:12:37 EDT
# mode: julia
    #scatter(vec(B[:,1]), vec(B[:,2]), title="$datacut")
# time: 2022-09-08 16:12:37 EDT
# mode: julia
    #Plot.save((;desc="scatter goal_index versus bestshift_bitsperspike",
# time: 2022-09-08 16:12:37 EDT
# mode: julia
    #           datacut))
# time: 2022-09-08 16:12:37 EDT
# mode: julia
    #histogram2d(vec(fwg[:bestshift_bitsperspike]), vec(gm), title="$datacut")
# time: 2022-09-08 16:12:37 EDT
# mode: julia
    #Plot.save((;desc="histogram2d goal_index versus bestshift_bitsperspike",
# time: 2022-09-08 16:12:37 EDT
# mode: julia
    #           datacut))
# time: 2022-09-08 16:12:37 EDT
# mode: julia
    #histogram(vec(fwg[:bestshift_bitsperspike]), title="$datacut")
# time: 2022-09-08 16:12:37 EDT
# mode: julia
    #Plot.save((;desc="histogram bestshift_bitsperspike", datacut))
# time: 2022-09-08 16:12:37 EDT
# mode: julia
    inds = sortperm(fwg[:bestshift_bitsperspike][:,1])
# time: 2022-09-08 16:12:37 EDT
# mode: julia
    bps  = fwg[:bitsperspike][inds, :]
# time: 2022-09-08 16:12:37 EDT
# mode: julia
    #X = hcat([Utils.norm_extrema(b) for b in eachrow(bps)]...)'
# time: 2022-09-08 16:12:37 EDT
# mode: julia
    #heatmap(sh, collect(1:size(X,1)), X, title="$datacut", size=(400,1000))
# time: 2022-09-08 16:12:37 EDT
# mode: julia
    #vline!([0],c=:black,linestyle=:dash, linewidth=2)
# time: 2022-09-08 16:12:37 EDT
# mode: julia
    #Plot.save((;desc="snake plot, norm 01", datacut))
# time: 2022-09-08 16:12:37 EDT
# mode: julia
    #X = hcat([Utils.norm_percent(b,0.5) for b in eachrow(bps)]...)'
# time: 2022-09-08 16:12:37 EDT
# mode: julia
    #heatmap(sh, collect(1:size(X,1)), X, title="$datacut", size=(400,1000))
# time: 2022-09-08 16:12:37 EDT
# mode: julia
    #vline!([0],c=:black,linestyle=:dash, linewidth=2)
# time: 2022-09-08 16:12:37 EDT
# mode: julia
    #Plot.save((;desc="snake plot, norm percent", datacut))
# time: 2022-09-08 16:12:37 EDT
# mode: julia
end
# time: 2022-09-08 16:12:37 EDT
# mode: julia
end

include(scriptsdir("timeshift","Tau-relation-to-goal.jl"))
future = sum(shiftdist .> 0 .&& abs.(shiftdist) .< 1)
past = sum(shiftdist .< 0 .&& abs.(shiftdist) .< 1)
k1,k2=:cue,:memory
@info "keys" k1 k2
key1,key2 = keyz[k1],keyz[k2]
f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]));
u = intersect(f1.dims[1], f2.dims[1])
vals = f1[unit=At(u)][:bestshift_bitsperspike],
           f2[unit=At(u)][:bestshift_bitsperspike]
edges = LinRange(minimum(f1[:shift]), maximum(f1[:shift]), bins+1)
xlim = (edges[begin], edges[end])
shiftdist = (vals[2] - vals[1])[:,1]
plot_histdist = 
        histogram(shiftdist; bins=edges,  xlim, normalize=:probability, alpha=0.5, label="", xlabel="shift", ylabel="pdf")
vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
dist = fit(Histogram, shiftdist, edges)
dist = StatsBase.normalize(dist, mode=:probability)
plot_histdist_cdf = 
        bar(dist.edges, cumsum(dist.weights); xlim, alpha=0.5, label="", xlabel="shift", ylabel="cdf")
vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
future = sum(shiftdist .> 0 .&& abs.(shiftdist) .< 1)
past = sum(shiftdist .< 0 .&& abs.(shiftdist) .< 1)
past

