using Plots

hatraj  = replace(beh.hatraj, missing=>"")
uHatraj = unique(skipmissing(beh.hatraj))

# CUEMEM frequency
H = []
for traj in uHatraj
    inds =  hatraj .== traj
    push!(H, Plots.histogram(beh[inds, :cuemem], bins=-1.5:0.5:1.5, title=traj, xticks=[-1,0,1]))
end
Plots.plot(H..., size=(1400,1400))

# CORRECT FREQUENCY
H = []
for traj in uHatraj
    inds =  hatraj .== traj
    push!(H, Plots.histogram(beh[inds, :correct], bins=-1.5:0.5:1.5, title=traj, xticks=[-1,0,1]))
end
Plots.plot(H..., size=(1400,1400))


