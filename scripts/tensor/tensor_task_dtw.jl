
using GoalFetchAnalysis
using Munge.tensor
using Munge.spiking
using Table
using Missings
using Plot.task

@time spikes, beh, ripples, cells = Load.load("RY16", 36);
tsk = Load.load_task("RY16",36);

#beh.G = collect(eachrow([beh.startWell beh.stopWell]))


#unique(beh.G)
unique(beh.subblock)
unique(beh.traj)
dims = ["startWell","stopWell","traj"]
values = ["x","y"]

X = torate(spikes, beh)
X = rate_todataframe(X, (beh,"time",[values...,dims...],))
X = Table.group.equalize(X,  setdiff(dims,["traj"]), :traj, thresh=8)

@time out = tensorize(X, dims, values)

using TensorToolbox
out = tenmat(out, 3)
mis = vec((!).(any(ismissing.(out), dims=1)))
out = out[:,findall(mis)]
out = disallowmissing(out)
out = [Matrix(o') for o in out]

"""
                                                  
,-.-.         |o              |                  o
| | |,---.,---|.,---.,---.    |--- ,---.,---.    .
| | ||---'|   ||,---||   |    |    |    ,---|    |
` ' '`---'`---'``---^`   '    `---'`    `---^    |
                                             `---'
"""

# APPROACH 1

# NOTE because it's hard to exactly grab an ideal trajectory, it could
# be an interesting idea to project all trajectories onto the line
# spanning from the start well to the stop well

using DynamicAxisWarping
function dtwmedian()
    col = 1
    I1 = Matrix{Vector}(undef, length(row), length(row))
    I2 = Matrix{Vector}(undef, length(row), length(row))
    C = Matrix{Vector}(undef, length(row), length(row))

    # Basically register each pattern to other patterns of the same type
    for (i,j) in Iterators.product(1:length(row), 1:length(row))
        cost, i1, i2 = dtw(row[i], row[j])
        I1[i,j] = i1
        I2[i,j] = i2
        C[i,j] = cost
    end

    # Compute the median pattern, by interpolating them all to the same size
    # and then taking the median
    for (i,j) in Iterators.product(1:length(row), 1:length(row))
    end

    # Perform dtw of all patterns to the median
    
    # And the cost-weighted mean
    
end

# APPROACH 2 : Included in DynamicAxisWarping.jl, Barycentric Averaging
r=1
P = []
for i in 1:size(out,2)
    p = plot()
    [plot!(eachrow(m)...;linestyle=:dash,c=:black) for m in reshape.(out[:,i], 2,:)]
    p
    Q = dba(vec.([c[1] for c in collect((eachrow(out[:,i])))]), DTW(r;transportcost=1.5))
    plot!(eachrow(reshape(Q[1],2,:))...)
    push!(P,p)
end
plot!(P[1], title="Barycentric")
plot!(P[2], title="Averaging\n=$r")
plot!(P[4], title="(definitely\nflawed)")
plot(P..., legend=:none, xticks=:none,yticks=:none)

# APPROACH 3 : Included in DynamicAxisWarping.jl, Soft Barycentric Averaging
bc = pyimport("tslearn.barycenters")
γ = 1
P = []
for i in 1:size(out,2)

    p=plotboundary(tsk)
    data = reshape.(out[:,i], 2,:)
    [plot!(eachrow(m)...;linestyle=:dash,c=:black) for m in data]
    p

    data = [Matrix(d') for d in data]
    dba_vecs = vec.([c[1] for c in collect((eachrow(out[:,i])))])
    #Q = dba(dba_vecs, SoftDTW(;γ,radius=r,transportcost=1.5))
    @time Q = bc.softdtw_barycenter(data, γ)
    #Q = reshape(Q,2,:)
    Q = hcat([rollmedian(q,2) for q in eachcol(Q)]...)
    plot!(eachcol(Q)...)
    push!(P,p)

end
plot!(P[1], title="Barycentric")
plot!(P[2], title="Averaging\n=$r")
plot!(P[4], title="(definitely\nflawed)")
plot(P..., legend=:none, xticks=:none,yticks=:none)


# APPROACH 4 : Literal MEDIAN
γ = 1
P = []
for i in 1:size(out,2)

    p=plot()
    data = reshape.(out[:,i], 2,:)
    maxdatalen = maximum([size(d,2) for  d in data])
    [plot!(eachrow(m)...;linestyle=:dash,c=:black) for m in data]
    p

    data = [Matrix(d') for d in data]

    # Interpolate
    for i in 1:length(data)
        D = collect.(eachcol(data[i]))
        D = interpolate.(D, [BSpline(Linear())])
        D = hcat([d[LinRange(1,length(d), maxdatalen)]
                  for d in D]...)
        data[i]=D
    end

    Q = median(cat(data...;dims=3),dims=3)[:,:]
    
    
    Q = hcat([rollmedian(q,2) for q in eachcol(Q)]...)
    plot!(eachcol(Q)...)
    push!(P,p)
end
plot!(P[1], title="Barycentric")
plot!(P[2], title="Averaging\n=$r")
plot!(P[4], title="(definitely\nflawed)")
plot(P..., legend=:none, xticks=:none,yticks=:none)
