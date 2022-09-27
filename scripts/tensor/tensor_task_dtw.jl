
using GoalFetchAnalysis
using Munge.tensor
using Munge.spiking

@time spikes, beh, ripples, cells = Load.load("RY16", 36);

#beh.G = collect(eachrow([beh.startWell beh.stopWell]))


unique(beh.G)
unique(beh.subblock)
unique(beh.traj)
dims = ["startWell","stopWell","subblock","traj"]
values = ["x","y"]

X = torate(spikes, beh)
X = rate_todataframe(X, (beh,"time",[values...,dims...],))

@time out = tensorize(X, dims, values)



