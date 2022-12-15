using GoalFetchAnalysis
using Munge.tensor

@time spikes, beh, ripples, cells = Load.load("RY16", 36);
dims = [:epoch, :traj]
vars = :x

# Quantilize works?
X = beh[rand(1:size(beh,1),1000),:]
x = copy(x)
quantilize(X, Dict{Symbol,Int}(:startWell=>3))
scatter(x.startWell, X.startWell)

# Relativize works?
X = beh[rand(1:size(beh,1),1000),:]
x = copy(x)
X = relativize(X, dims, :time; keeporig=true)
scatter(X.origtime, X.time)

# Tensorize works?
X = beh
@time xy = tensor_continuous(X, [:traj], [:y,:x])
@time x  = tensor_continuous(X, [:traj],  [:x])
@time y  = tensor_continuous(X, [:traj],  [:y])

# ============
# Trial x Goal
# ============

dims = [:trial, :stopwell]
val_beh  = [:x, :y]
val_rate = :rate
# Get the behavior tensor and equalize traj
B = tensor_continuous(beh, dims, val_beh)
B = equalize(B, :traj)
# Dynamic warp
templates       = tmedian(B, :traj)
# expectation-maximize templates
for i in 1:10
    B̂, cost, i1, i2 = tdtw(templates, B)
    templates       = tmedian(B, :traj)
end
# Visualize


S = tensor.tensor_pointproc(spikes, dims, val_rate)
S = apply_tdtw(i2, S, B)
rate_templates = tmedian(S, :traj)
# Visualize
