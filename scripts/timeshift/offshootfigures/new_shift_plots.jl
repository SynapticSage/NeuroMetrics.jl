using GoalFetchAnalysis
using Serialization
using Plots
using Timeshift.shiftmetrics
using Field.metrics
using Timeshift
using DimensionalData
using Statistics, NaNStatistics, Bootstrap
import Plot
using StatsBase
using ProgressMeter
using DataFramesMeta

F,I,shifts = deserialize(datadir("fixed_shifts_smallerbounds_higherres_2.serial"));
keyz = Dict(key.datacut=>key for key in keys(F))

parfolder="highres2"
theme(:bright)
function prep(f)    
    push_dims!(f)
    pop_metric!(f, :unit)
    push_metric!(f, metrics.bitsperspike)
    push_metric!(f, metrics.totalcount)
    push_metric!(f, metrics.maxrate)
    push_shiftmetric!(f, best_tau!; metric=:bitsperspike)
    push_metric!(f, metrics.bitsperspikeold)
    push_shiftmetric!(f, best_tau!; metric=:bitsperspikeold)
    f = f[vec(all(f[:totalcount] .> 50, dims=2) .&&
              any(f[:bitsperspike] .> 0.5, dims=2)) ,:]
end
  
function doheat(f,key)
    inds = sortperm(f[:bestshift_bitsperspike][:,1])
    bps  = f[:bitsperspike][inds, :]
    XX = hcat([Utils.norm_extrema(b) for b in eachrow(bps)]...)'
    heatmap(shifts, collect(1:size(XX,1)), XX, title="$(key.datacut)", size=2 .* (600,800))
end


Plot.setfolder(parfolder,"hist_diff_correct")
bins = 25
for (k1, k2) in zip((:nontask, :cue), (:task,:memory))
@info "keys" k1 k2
key1,key2 = keyz[k1],keyz[k2]
    f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]));
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
         
     future = sum(shiftdist .> 0 .&& abs.(shiftdist) .< 1)
     past = sum(shiftdist .< 0 .&& abs.(shiftdist) .< 1)
                     
     plot(plot_histdist, plot_histdist_cdf,
          bar(["past","future"],[past,future]./(past+future)),
           layout=[grid(2,1), grid(1,1)])
           
   #scatter(vec.(vals)...)
   #plot!(-2:2,-2:2,c=:white)
   Plot.save((;k1,k2))
end

func = mean
Plot.setfolder(parfolder,"fixed_$(string(Symbol(func)))maps3")
# time: 2022-09-08 18:05:44 EDT
# mode: julia
# Visualize the heatmap
# time: 2022-09-08 18:05:44 EDT
# mode: julia
P=[]
# time: 2022-09-08 18:05:45 EDT
# mode: julia
for (k1, k2) in zip((:nontask, :cue, :mem_correct, :cue_correct, :correct), (:task,:memory, :mem_error, :cue_error, :error))
    @info "keys" k1 k2
    key1,key2 = keyz[k1],keyz[k2]
    f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]))
    u = intersect(f1.dims[1], f2.dims[1])
    f1,f2 = f1[unit=At(u)], f2[unit=At(u)]

    inds = sortperm(f1[:bestshift_bitsperspike][:,1])
    XX1 =bps  = f1[:bitsperspike][inds, :]
    V1 = vec(func(XX1,dims=1))
    bs1 = bootstrap(func, collect(eachrow(XX1)), BalancedSampling(1000))
    straps1 = hcat(straps(bs1)...)

    inds = sortperm(f2[:bestshift_bitsperspike][:,1])
    XX2 = bps  = f2[:bitsperspike][inds, :]
    V2 = vec(func(XX2,dims=1))
    bs2 = bootstrap(func, collect(eachrow(XX2)), BasicSampling(1000))
    straps2 = hcat(straps(bs1)...)

    # I can do this because I intersected neurons in both
    bs21 = bootstrap(func, collect(eachrow(XX2-XX1)), BasicSampling(1000))
    ci21 = confint(bs21, BasicConfInt(0.95))
    straps21 = [x-y for (x,y) in Iterators.product(eachrow(straps1),
                                                   eachrow(straps2))]
    shifts = collect(f1.dims[2])
    # Manual jackknife
    M = Matrix(XX2-XX1)
    jacknives = [mean(M[setdiff(1:size(M,1),[i]),:],dims=1)
                 for i in 1:size(M,1)]
    plot(hcat([[y for y in x] for x in ci21]...)')
    
    p=plot(shifts, V2-V1, title="$k2 - $k1", size=0.75 .* (800,1200), fillrange=0,
          xlabel="shift", ylabel="Δ bits per spike")
    vline!([0],c=:black,linestyle=:dash, linewidth=2)
    plot!(shifts,hcat([vec(m) for m in jacknives]...),c=:black,legend=:none, linestyle=:dot, alpha=0.4,title="Jacknifed $(string(Symbol(func))) difference in\n bits per spike, $k2-$k1");
    hline!([0],c=:red,linestyle=:dash)
    
    if (k1,k2) == (:cue, :memory)
        shifts_local = Utils.in_range(shifts, [-.21, .61])
        shifts_distal = (!).(shifts_local)
        get_local(val)  = nanmean(val[shifts_local])
        get_distal(val) = nanmean(val[shifts_distal])
        L, D = get_local((V2-V1)),
               get_distal((V2-V1))
       stde(x) = nanstd(x,dims=1)/sqrt(size(x,1))
       J = vcat(jacknives...)
       Lse, Dse = mean(get_local(stde(J))), mean(get_distal(stde(J)))
       b=bar(["local", "nonlocal"],   [L,D], yerror=[Lse, Dse], label="")
       ##scatter!(["local", "nonlocal"], [L,D], err=[Lse, Dse],
                #markerstrokecolor = :transparent,
                #linestyle=:dash, c=:black, markersize=10, label="", lw=4)
       plot!(["local", "nonlocal"], [L,D], err=[Lse, Dse],
                l=nothing, m=(4,stroke(:red, 1.5)), label="")

       
       plot(p, b, size=(1000,1000))
       Plot.save((;desc="barsum", k1, k2))
    end
    push!(P,p)
    p
    Plot.save((;xlim="full", k1, k2))
end
# time: 2022-09-08 18:05:51 EDT
# mode: julia
plot(P...)
# time: 2022-09-08 18:05:51 EDT
# mode: julia
Plot.save("together")

func = mean
nanfunc = eval(Symbol("nan"*string(Symbol(func))))
Plot.setfolder(parfolder,"fixed_$(string(Symbol(func)))maps_centered")
# time: 2022-09-08 18:05:44 EDT
# mode: julia
# Visualize the heatmap
# time: 2022-09-08 18:05:44 EDT
# mode: julia
P=[]
# time: 2022-09-08 18:05:45 EDT
# mode: julia
for (k1, k2) in zip((:nontask, :cue, :mem_error, :cue_error, :error), (:task,:memory, :mem_correct, :cue_correct, :correct))
    @info "keys" k1 k2
    key1,key2 = keyz[k1],keyz[k2]
    f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]))
    u = intersect(f1.dims[1], f2.dims[1])
    f1,f2 = f1[unit=At(u)], f2[unit=At(u)]

    inds = sortperm(f1[:bestshift_bitsperspike][:,1])
    XX1 = bps  = f1[:bitsperspike][inds, :]
    XX1 = bps = XX1 .- Matrix(func(XX1, dims=2))
    V1 = vec(func(XX1,dims=1))
    #bs1 = bootstrap(func, collect(eachrow(XX1)), BalancedSampling(1000))
    #straps1 = hcat(straps(bs1)...)

    inds = sortperm(f2[:bestshift_bitsperspike][:,1])
    XX2 = bps  = f2[:bitsperspike][inds, :]
    XX2 = XX2 .- Matrix(func(XX2, dims=2))
    V2 = vec(func(XX2,dims=1))
    #bs2 = bootstrap(func, collect(eachrow(XX2)), BasicSampling(1000))
    #straps2 = hcat(straps(bs1)...)

    # I can do this because I intersected neurons in both
    #bs21 = bootstrap(func, collect(eachrow(XX2-XX1)), BasicSampling(1000))
    #ci21 = confint(bs21, BasicConfInt(0.95))
    #straps21 = [x-y for (x,y) in Iterators.product(eachrow(straps1),
                                                   #eachrow(straps2))]
    shifts = collect(f1.dims[2])
    # Manual jackknife
    M = Matrix(XX2-XX1)
    jacknives = [func(M[setdiff(1:size(M,1),[i]),:],dims=1)
                 for i in 1:size(M,1)]
    plot(hcat([[y for y in x] for x in ci21]...)')
    
    p=plot(shifts, V2-V1, title="$k2 - $k1", size=0.75 .* (800,1200), fillrange=0,
          xlabel="shift", ylabel="Δ bits per spike")
    vline!([0],c=:black,linestyle=:dash, linewidth=2)
    plot!(shifts,hcat([vec(m) for m in jacknives]...),c=:black,legend=:none, linestyle=:dot, alpha=0.4,title="Jacknifed $(string(Symbol(func))) difference in\n bits per spike, $k2-$k1");
    hline!([0],c=:red,linestyle=:dash)
    
    if (k1,k2) == (:cue, :memory) && func == mean
        shifts_local = Utils.in_range(shifts, [-.21, .61])
        shifts_distal = (!).(shifts_local)
        get_local(val)  = nanfunc(val[shifts_local])
        get_distal(val) = nanfunc(val[shifts_distal])
        L, D = get_local((V2-V1)),
               get_distal((V2-V1))
       stde(x) = nanstd(x,dims=1)/sqrt(size(x,1))
       J = vcat(jacknives...)
       Lse, Dse = func(get_local(stde(J))), func(get_distal(stde(J)))
       b=bar(["local", "nonlocal"],   [L,D], yerror=[Lse, Dse], label="")
       ##scatter!(["local", "nonlocal"], [L,D], err=[Lse, Dse],
                #markerstrokecolor = :transparent,
                #linestyle=:dash, c=:black, markersize=10, label="", lw=4)
       plot!(["local", "nonlocal"], [L,D], err=[Lse, Dse],
                l=nothing, m=(4,stroke(:red, 1.5)), label="")

       
       plot(p, b, size=(1000,1000))
       Plot.save((;desc="barsum", k1, k2))
    end
    push!(P,p)
    p
    Plot.save((;xlim="full", k1, k2))
end
# time: 2022-09-08 18:05:51 EDT
# mode: julia
plot(P...)
# time: 2022-09-08 18:05:51 EDT
# mode: julia
Plot.save("together")

# Seems to uncover an effect in the mean bits per spike,
# with current being -0.2s to 0.4s


Plot.setfolder(parfolder,"fixed_heatmaps3")
for key in keys(F)
    f = prep(matrixform(F[key]))
    #push_metric!(f, metrics.bitsperspike)
    #push_shiftmetric!(f, best_tau!; metric=:bitsperspike)
    shifts = vec(f.dims[2].val.data)
    inds = sortperm(f[:bestshift_bitsperspike][:,1])
    bps  = f[:bitsperspike][inds, :]
    XX = hcat([Utils.norm_extrema(b) for b in eachrow(bps)]...)'
    heatmap(shifts, collect(1:size(XX,1)), XX, 
    title="$(key.datacut)", size=(600,1200), c=:hot, clim=(0.3,1))
    vline!([0],c=:black,linestyle=:dash, linewidth=2)
    Plot.save((;xlim="full", key...))
end


function field_grad(f)
    np = pyimport("numpy")
    ff = [ff.rate for ff in f]
    cell = cat(ff[2,:]...,dims=3)
    np.gradient(cell)
end

function viz_gradient(nabla, scale=1, ck=false, scalek=false)

    xx,yy,zz,ii,jj,kk = ([] for _ in 1:6)
    for (xyz, i, j, k) in zip(CartesianIndices(nabla[1]), nabla...)
        x,y,z = xyz.I
        sc = scale isa Vector ? scale[z] : scale
        push!(xx,x)
        push!(yy,y)
        push!(zz,z)
        push!(ii,i*scale)
        push!(jj,j*scale)
        push!(kk,k*scale)
    end
    plt = pyimport("matplotlib.pyplot")
    #fig,axs = plt[:subplots](Int(ceil(sqrt(maximum(zz)))), Int(ceil(sqrt(maximum(zz)))))
    
    @gif for zi in 1:160
        x,y,z,i,j, k =getindex.([xx,yy,zz,ii,jj,kk], [zz.==zi])
        #axs[zi][:quiver](x, y, i, j)
        arrw = (i, j)
        arrw = scalek ? arrw .* abs.(k) : arrw
        quiver(x,y, quiver= 0.05 .* arrw , c=:black)
    end

end


using CausalityTools
@time spikes, beh, ripples, cells = Load.load("RY16", 36);
Rca1, Rpfc = (Munge.spiking.torate(@subset(spikes,:area .== ar), beh) for ar in ("CA1","PFC"))

P=[]
prog = Progress(size(Rca1,2)*size(Rpfc,2); desc="transfer ent")
sets = collect(Iterators.product(1:size(Rca1,2), 1:size(Rpfc,2)))
P = Matrix{Vector}(undef,size(sets))
Q = Matrix{Vector}(undef,size(sets))
D = falses(size(sets))
using Threads

@showprogress for s in sets
    i,j = s
    if D[i,j] == false
        p = Threads.@spawn predictive_asymmetry(Rca1[:,i], Rpfc[:,j], VisitationFrequency(RectangularBinning(5)), 1:10)
        q = Threads.@spawn predictive_asymmetry(Rpfc[:,j], Rca1[:,i], VisitationFrequency(RectangularBinning(5)), 1:10)
        P[i,j], Q[i,j] = fetch(p), fetch(q)
        next!(prog)
        D[i,j]=true
    end
end

mkpath(datadir("transent"))
serialize(datadir("transent","alltimes.serial"),(;P,Q,D))

p=plot();[plot!(ppp,c=:black,linestyle=:dash,alpha=0.1) for ppp in pp]; p

# ==========================================

sets = collect(Iterators.product(1:size(Rca1,2), 1:size(Rpfc,2)))
beh.cuemem = Int16.(beh.cuemem)
pa    = Dict{Tuple, Matrix{Vector}}()
done  = Dict{Tuple, Matrix{Bool}}()
savefile = datadir("transent","condtimes.serial")
checkpoint = 100

# TODO create a traj 3-bin (still at well, 1st half move, 2nd half move)
@showprogress "splits" for (cuemem, corr) in Iterators.product(unique(beh.cuemem), unique(beh.correct))
    # Select data
    ind = beh.cuemem .== cuemem .&& beh.correct .== corr
    if sum(ind) == 0; continue; end
    rca1, rpfc = Rca1[ind,:], Rpfc[ind,:]
    # Initialize data trackers
    pa[(cuemem, corr, "pfc->ca1")] = Matrix{Vector}(undef,size(sets))
    pa[(cuemem, corr, "ca1->pfc")] = Matrix{Vector}(undef,size(sets))
    done[(cuemem, corr)] = falses(size(sets))
    d = done[(cuemem, corr)]
    prog = Progress(size(Rca1,2)*size(Rpfc,2); desc="sets")
    for (count, s) in enumerate(sets)
        i,j = s
        if d[i,j] == false
            x, y = vec.((rca1[:,i], rpfc[:,j]))
            #xc, yc = CuArray.(vec.((rca1[:,i], rpfc[:,j])))
            p = Threads.@spawn predictive_asymmetry(x, y, 
                          VisitationFrequency(RectangularBinning(5)), 1:10)
            #@time p = predictive_asymmetry(xc, yc, 
                                           #VisitationFrequency(RectangularBinning(5)), CuArray(1:10))
            q = Threads.@spawn predictive_asymmetry(y, x, 
                          VisitationFrequency(RectangularBinning(5)), 1:10)
            pa[(cuemem, corr, "pfc->ca1")][i,j], 
            pa[(cuemem, corr, "pfc->ca1")][i,j] = Vector(fetch(p)), Vector(fetch(q))
            next!(prog)
            d[i,j]=true
        end
        mod(count, checkpoint) == 0 ? serialize(condtimes, (;pa, done)) : nothing
    end
    serialize(condtimes, (;pa, done))
end

# ==========================================

using UMAP
@time spikes, beh, ripples, cells = Load.load("RY16", 36);
Rca1, Rpfc = (Munge.spiking.torate(@subset(spikes,:area .== ar), beh) for ar in ("CA1","PFC"))

ca1_embedding = umap((Rca1')[:,1:1000], 3)
pfc_embedding = umap(Rpfc, 3)
