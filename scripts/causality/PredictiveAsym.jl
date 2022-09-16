
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

using CausalityTools
@time spikes, beh, ripples, cells = Load.load("RY16", 36);
Rca1, Rpfc = (Munge.spiking.torate(@subset(spikes,:area .== ar), beh) for ar in ("CA1","PFC"))

P=[]
prog = Progress(size(Rca1,2)*size(Rpfc,2); desc="transfer ent")
sets = collect(Iterators.product(1:size(Rca1,2), 1:size(Rpfc,2)))
P = Matrix{Vector}(undef,size(sets))
Q = Matrix{Vector}(undef,size(sets))
D = falses(size(sets))
Dq = reshape([isassigned(Q,p) for p in eachindex(Q)],size(Q))
using Threads

@showprogress for s in sets
    i,j = s
    if Dq[i,j] == false
        #p = Threads.@spawn predictive_asymmetry(Rca1[:,i], Rpfc[:,j], VisitationFrequency(RectangularBinning(5)), 1:10)
        q = predictive_asymmetry(Rpfc[:,j], Rca1[:,i], VisitationFrequency(RectangularBinning(5)), 1:10)
        #P[i,j], Q[i,j] = fetch(p), fetch(q)
        next!(prog)
        Dq[i,j]=true
    end
end

mkpath(datadir("transent"))
serialize(datadir("transent","alltimes.serial"),(;P,Q,D))

qq = Q[[isassigned(Q,p) for p in eachindex(Q)]]
qq = hcat(qq...)'

p=plot();[plot!(ppp,c=:black,linestyle=:dash,alpha=0.1) for ppp in pp]; p
p=plot();[plot!(qqq,c=:black,linestyle=:dash,alpha=0.1) for qqq in qq]; p

plot!(mean(qq,dims=1)', m=:circle, fillrange=0, label="pfc->ca1")
hline!([0])


plot(mean(abs.(pp),dims=1)', m=:circle)
plot!(mean(abs.(qq), dims=1)', m=:circle)
hline!([0])

@userplot PlotMeanCause
@recipe function plotmeancause(plt::PlotMeanCause; transform=identity)
    P = plt.args[1]
    pp = P[[isassigned(P,p) for p in eachindex(P)]]
    pp = hcat(pp...)'
    m --> :circle
    fillrange --> 0
    label --> "thing a -> thing b"
    #@infiltrate
    mean(transform(pp),dims=1)'
end

# ==========================================

sets = collect(Iterators.product(1:size(Rca1,2), 1:size(Rpfc,2)))
beh.cuemem = Int16.(beh.cuemem)
pa    = Dict{Tuple, Matrix{Vector}}()
done  = Dict{Tuple, Matrix{Bool}}()
savefile = datadir("transent","condtimes.serial")
checkpoint = 100
for (cuemem, corr) in Iterators.product(unique(beh.cuemem), unique(beh.correct))
    done[(cuemem, corr)] = falses(size(sets))

# TODO create a traj 3-bin (still at well, 1st half move, 2nd half move)
@showprogress "splits" for (cuemem, corr) in Iterators.product(unique(beh.cuemem), unique(beh.correct))
    @info (cuemem, corr)
    if (cuemem == -1 && corr == -1) ||
        isnan(cuemem) || isnan(corr)
        continue
    end
    # Select data
    ind = beh.cuemem .== cuemem .&& beh.correct .== corr
    if sum(ind) == 0; continue; end
    rca1, rpfc = Rca1[ind,:], Rpfc[ind,:]
    # Initialize data trackers
    pa[(cuemem, corr, "pfc->ca1")] = Matrix{Vector}(undef,size(sets))
    pa[(cuemem, corr, "ca1->pfc")] = Matrix{Vector}(undef,size(sets))

    d = done[(cuemem, corr)]
    prog = Progress(size(Rca1,2)*size(Rpfc,2); desc="sets")
    for (count, s) in enumerate(sets)
        i,j = s
        if d[i,j] == false
            x, y = vec.((rca1[:,i], rpfc[:,j]))
            if all(x.==0) || all(y.==0) 
                p, q = zeros(10), zeros(10)
            #xc, yc = CuArray.(vec.((rca1[:,i], rpfc[:,j])))
            else
                x, y = zscore(x), zscore(y)
                p = Threads.@spawn predictive_asymmetry(x, y, VisitationFrequency(RectangularBinning(5)), 1:10)
                #@time p = predictive_asymmetry(xc, yc, 
                                               #VisitationFrequency(RectangularBinning(5)), CuArray(1:10))
                q = Threads.@spawn predictive_asymmetry(y, x, VisitationFrequency(RectangularBinning(5)), 1:10)
            end
            try
                pa[(cuemem, corr, "ca1->pfc")][i,j], 
                pa[(cuemem, corr, "pfc->ca1")][i,j] = Vector(fetch(p)), Vector(fetch(q))
            catch
                @infiltrate
            end
            next!(prog)
            d[i,j]=true
        end
        mod(count, checkpoint) == 0 ? serialize(savefile, (;pa, done)) : nothing
    end
    serialize(savefile, (;pa, done))
end

# ==========================================

label = "pfc->ca1"
plotmeancause(pa[(1, 1, label)])

