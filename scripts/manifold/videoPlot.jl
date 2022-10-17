
using Images
using Makie, GLMakie
video = Load.video.load_videocollection("RY16", 36)
times = beh.time[1:size(em,1)]

t = 10000
times = times[t]

prog = Progress(length(times);desc="Times")
anim = @animate for (t,time) in enumerate(times)
    @time p1 = plot(eachcol(em)...; alpha=0.01, c=:gray)  
    scatter!(eachrow(em[t,:])...; c=:red, markersize=4, legend=:bottomleft, label="current")
    @time p2 = plot(video[times[1]], legend=:none)
    plot(p1,p2, size=(1000,600))
    next!(prog)
end

