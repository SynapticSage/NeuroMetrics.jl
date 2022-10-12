
using DynamicAxisWarping
using GoalFetchAnalysis
using Munge.tensor
using Munge.spiking
using Table
using Missings
using Plot.task
using Infiltrator
using RecipesBase

@time spikes, beh, ripples, cells = Load.load("RY16", 36);
tsk = Load.load_task("RY16",36);

#unique(beh.G)
unique(beh.subblock)
unique(beh.traj)
dims = ["startWell","stopWell","traj"]
values = ["x","y", "time"]
n = length(values)-1

X = torate(spikes, beh)
X = rate_todataframe(X, (beh,"time",[values...,dims...],))
X = Table.group.equalize(X,  setdiff(dims,["traj"]), :traj, thresh=12)

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

begin
function dtwmedian()
    col = 1
    I1 = Matrix{Vector}(undef, length(row), length(row))
    I2 = Matrix{Vector}(undef, length(row), length(row))
    C  = Matrix{Vector}(undef, length(row), length(row))

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
end

# APPROACH 2 : Included in DynamicAxisWarping.jl, Barycentric Averaging
begin
r=1
P = []
for i in 1:size(out,2)
    p = plot()
    [plot!(eachrow(m)...;linestyle=:dash,c=:black) for m in out[1:n,:]]
    p
    Q = dba(vec.([c[1] for c in collect((eachrow(out[:,i])))]), 
            DTW(r;transportcost=1.5))
    plot!(eachrow(reshape(Q[1],n,:))...)
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
    data = out[:,i][1:n,:]
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
end


# APPROACH 4 : Literal MEDIAN
begin
P = []
templates = []
for i in 1:size(out,2)

    p=plot()
    data = out[:,i][1:n,:]
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

    template = median(cat(data...;dims=3),dims=3)[:,:]
    template = hcat([rollmedian(q,2) for q in eachcol(template)]...)
    push!(templates, template)
    plot!(eachcol(template)...)

    push!(P,p)
end
plot!(P[1], title="Barycentric")
plot!(P[2], title="Averaging\n=$r")
plot!(P[4], title="(definitely\nflawed)")
plot(P..., legend=:none, xticks=:none,yticks=:none)
end

# APPLY TEMPLATES :: DO DTW
templag = []
correction = []
for i in 1:size(out,2)

    data = out[:,i]
    inputdata = getindex.(out[:,i], [1:n], [:])
    dtwres = dtw.([templates[i][:,1:n]'], inputdata)
    seq1 = [a[2] for a in dtwres]
    seq2 = [a[3] for a in dtwres]
    push!(correction, seq2)
    push!(templag, seq1)

    # Apply warps
    datawarp = [w[:,s2]   for (w,s2) in zip(inputdata, seq2)]
    timewarp = [w[end,s2] for (w,s2) in zip(data, seq2)]
    
    template_time = templates[i][:,end]
    template_time .-= minimum(template_time)

    # Get time effects
    warptable = [DataFrame([template_time[s1] d[end, s2] s1 s2],["time_template","time","s1","s2"])
                 for (s1,s2,d) in zip(seq1, seq2, data)]

    p1=plot()
    [plot!(eachrow(m)...;linestyle=:dash,c=:black) for m in inputdata]
    plot!(eachrow(templates[i]')...,linewidth=3,linestyle=:solid)
    p2=plot()
    [plot!(eachrow(m)...;linestyle=:dash,c=:black) for m in datawarp]
    plot!(eachrow(templates[i]')...,linewidth=3,linestyle=:solid)
    p3=plot()
    [plot!(s,fillrange=0,legend=:none) for s in seq2]
    plot(p1,p2,p3, layout=[grid(1,2); grid(1,1)])

end

@userplot MatchPlot2
znorm2(x) = (x = x.- mean(x,dims=2); x ./= std(x,dims=2))
@recipe function f(h::MatchPlot2; 
                   transportcost=1, separation=0.5, ds=1, postprocess=nothing)

    x, y, D, i1, i2 = DynamicAxisWarping.handleargs(h;
                                                    transportcost=transportcost,
                                                    postprocess=postprocess)
    x,y = znorm2.((x,y))
    s1 = x .- separation
    s2 = y .+ separation

    @series begin
        @infiltrate
        (collect(eachrow(s1))...,)
    end
    @series begin
        (collect(eachrow(s2))...,)
    end
    @series begin
        primary := false
        linecolor --> :black
        seriesalpha --> 0.2
        s1, s2 = s1[:,i1][:,1:ds:end], s2[:,i2][:,1:ds:end]
        # Concatonate along 3rd dim
        s3 = cat(s1,s2; dims=3)
        s3 = eachrow.(eachslice(s3, dims=2))
        [(rows_of_slices...,) for rows_of_slices in s3]
    end
end

P = []
@showprogress for i in 1:size(out,2)
    data = reshape.(out[:,i], 2,:)
    M = [matchplot2(templates[i]',data[j]; ds=5)
         for j in 1:length(data)]
    push!(P,plot(M...))
end


using Munge.dynamic
examples = get_groupedexamples(spikes, beh)
templates = get_templates(examples)
dtwtab = get_dtwtable(examples, templates)

apply_warps(X, dtwtab)
