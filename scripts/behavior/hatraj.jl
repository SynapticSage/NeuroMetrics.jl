using Plots
import Plot

using GoalFetchAnalysis
using DI

animal, day = "RY16", 36
beh = load_behavior(animal, day)

Plot.setparentfolder()
Plot.setfolder("behavior", "hatraj")
Plot.setappend("_$animal.$day")

behsort = sort(beh, :hatraj)
hatraj  = replace(behsort.hatraj, missing=>"")
uHatraj = unique(skipmissing(behsort.hatraj))

# CUEMEM frequency
begin
    H = []
    for traj in uHatraj
        inds =  hatraj .== traj
        push!(H, Plots.histogram(behsort[inds, :cuemem], bins=-1.5:0.5:1.5, title=traj, xticks=[-1,0,1],label=""))
    end
    Plots.plot(H..., size=(1400,1400))
    Plot.save("cuemem frequency")
end

# CORRECT FREQUENCY
begin
    H = []
    for traj in uHatraj
        inds =  hatraj .== traj
        push!(H, Plots.histogram(behsort[inds, :correct], bins=-1.5:0.5:1.5, title=traj, xticks=[-1,0,1],label=""))
    end
    Plots.plot(H..., size=(1400,1400))
    Plot.save("correct frequency")
end

begin
    cm = sort(countmap(beh.hatraj))
    bar(cm; xticks=(collect(1:length(cm)) .- 0.5, collect(keys(cm))))

    b = subset(beh, :cuemem => c -> c .== 0, :correct => c-> c.==1)
    cm = sort(countmap(b.hatraj))
    b1=bar(cm; title="cue correct", xticks=(collect(1:length(cm)) .- 0.5, collect(keys(cm))))

    b = subset(beh, :cuemem => c -> c .== 1, :correct => c -> c .== 1)
    cm = sort(countmap(b.hatraj))
    b2=bar(cm; title="mem correct" ,xticks=(collect(1:length(cm)) .- 0.5, collect(keys(cm))))

    plot(b1,b2, layout=grid(2,1))

    Plot.save("By hatraj label, mem and cue correct")
end
