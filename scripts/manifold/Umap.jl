using DataStructures: OrderedDict

using GoalFetchAnalysis
import Utils.namedtup: ntopt_string
using  DataFramesMeta
using UMAP
@time spikes, beh, ripples, cells = Load.load("RY16", 36);
Rca1, Rpfc = (Munge.spiking.torate(@subset(spikes,:area .== ar), beh) for ar in ("CA1","PFC"))
nsamp = min(100_000, size(Rca1,1))

embedding, embedding2 = Dict(), Dict()
areas = (:ca1,:pfc)
embedding[:ca1] = umap((Rca1')[:,1:nsamp], 3)
embedding[:pfc] = umap((Rpfc')[:,1:nsamp], 3)
embedding2[:ca1] = umap((Rca1')[:,1:nsamp], 2)
embedding2[:pfc] = umap((Rpfc')[:,1:nsamp], 2)
embedding  = Dict(k=>v' for (k,v) in embedding)
embedding2 = Dict(k=>v' for (k,v) in embedding2)
inds = Dict()
for area in areas
    inds[area] = Utils.clean.inds_quantile_filter_dims(embedding[area], [0.02, 0.96])
end
serialize(datadir("manifold","ca1pfc_manifolds.serial"), (;embedding, embedding2, inds, animal="RY16",day=36))

# Stereoscopic views
for bfield in (:cuemem, :correct, :stopWell, :startWell, :trajType)
    cfull =  Utils.plot.plotcolor(beh[:,bfield],:vik)
    Plots.setfolder("manifold",string(bfield))
    for area in areas
        title="area=$area, bfield=$bfield"
        ia = inds[area]
        em = embedding[area][ia,:]
        em2 = embedding2[area][ia,:]
        c = cfull[1:nsamp][ia]

        @gif for i in 1:360; 
            plot(
                     scatter(eachcol(em)...; c=, camera=(mod(i+5,360),35), projection_type=:perspective),
                             scatter(eachcol(em)...; c=c[1:nsamp][ia], camera=(i,35), projection_type=:perspective);
                             title
                                        )
        end

        @gif for i in 1:360; 
            plot(
                     plot(eachcol(em)...; c=c[1:nsamp][ia], camera=(mod(i+5,360),35), projection_type=:perspective),
                             plot(eachcol(em)...; c=c[1:nsamp][ia], camera=(i,35), projection_type=:perspective);
                             title, alpha=0.05)
        end


        scatter(eachcol(em2)...; c, alpha=0.05, title)
        Plot.save("scatter_area=$area")
        plot(eachcol(em2)...; c, alpha=0.1, title)
        Plot.save("scatter_area=$area")

    end

end


#  ================
#  ISOLATED SPIKING
#  ================


# iso spikes
beh.index = 1:size(beh,1)
Utils.filtreg.register(beh,spikes,on="time",transfer=["index"])
lfp = Load.load_lfp("RY16", 36, tet=5);
lfp.time = lfp.time .- Load.min_time_records[1]
lfp = Munge.lfp.annotate_cycles(lfp)
Munge.spiking.isolated(spikes, lfp, include_samples=false)
spikes.super_isolated = spikes.isolated  .&& spikes.phase .> 3.14

I = combine(groupby(spikes,:index), :isolated=>mean, :super_isolated=>mean)[1:nsamp,:]
E = sqrt.(sum(diff(embedding[:ca1], dims=1).^2,dims=2))
E = [0; E]
I[!,:E] .= E

using GLM

using  StatsBase
using  StatsPlots
full = lm(@formula(E ~ isolated_mean), I)
#full = lm(@formula(E ~ super_isolated_mean), I)
adjr²(full)
@df I scatter(:isolated_mean, :E)

#lm(@formula(E ~ isolated_mean * isolated_mean^2 + isolated_mean^3), 
#   @subset(I, :isolated_mean .>0, :isolated_mean .< 1))
#lm(@formula(E ~ super_isolated_mean * super_isolated_mean^2 + super_isolated_mean^3), 
#   @subset(I, :super_isolated_mean .>0, :super_isolated_mean .< 1))

using Polynomials
@df @subset(I, :isolated_mean .>0, :isolated_mean .< 1) scatter(:isolated_mean, :E, alpha=0.01, aspect_ratio=0.01)
hline!([0], c=:black)

I.ibin = Utils.binning.digitize(I.isolated_mean, 10)
I.exist = I.isolated_mean .> 0
summ = combine(groupby(@subset(I, :isolated_mean .>=0, :isolated_mean .< 0.55), :ibin), 
        :E=>mean, 
        :E=>(x->std(x)/sqrt(length(x))) => :E_stde)

b=@df summ bar(:ibin, :E_mean, yerror=:E_stde, title="Manifold jump, iso  spikes",label="")
for i in 1:size(summ,1)
    #scatter!(b,[summ.ibin[i]], [summ.E_mean[i]], yerr=summ.E_stde[i],
             #lw=5, markerstrokecolor = :auto)
     val = [summ.E_mean[i]] .+ ([-0.5,0.5] .* summ.E_stde[i])
     @info val
     plot!([summ.ibin[i], summ.ibin[i]], val, c=:black, linewidth=5, linestyle=:solid, label="")
end
plot!()

summ = combine(groupby(@subset(I, :isolated_mean .>=0, :isolated_mean .< 0.55), :exist), 
        :E=>mean, 
        :E=>median, 
        :E=>(x->std(x)/sqrt(length(x))) => :E_stde)

b=@df summ bar(:exist, :E_mean, yerror=:E_stde, title="Manifold jump, iso  spikes",label="")
for i in 1:size(summ,1)
    #scatter!(b,[summ.ibin[i]], [summ.E_mean[i]], yerr=summ.E_stde[i],
             #lw=5, markerstrokecolor = :auto)
     val = [summ.exist[i]] .+ ([-0.5,0.5] .* summ.E_stde[i])
     @info val
     plot!([summ.ibin[i], summ.ibin[i]], val, c=:black, linewidth=5, linestyle=:solid, label="")
end
plot!()

#  ================
# CAUSALITY
#  ================

embedding_rows = Dict(k=>Dataset(v) for (k,v) in embedding);
est = VisitationFrequency(RectangularBinning(9))
params = (m = 15, τ = 1, horizon = 1:1000)

#params = (m = 100, τ = 4)
est = SymbolicPermutation(;params...)


cp_manifold_pa = Threads.@spawn CausalityTools.predictive_asymmetry(embedding_rows[:ca1],
                                                     embedding_rows[:pfc], est,
                                                     params.horizon);
pc_manifold_pa = Threads.@spawn CausalityTools.predictive_asymmetry(embedding_rows[:pfc],
                                                     embedding_rows[:ca1], est,
                                                     params.horizon);

serialize(datadir("manifold","manifold_pa_cause$(ntopt_string((;params...)))"),
                  (;params..., est, cp_manifold_pa, pc_manifold_pa)
                 )

task_vars = replace([beh.cuemem beh.correct],NaN=>-1)[1:nsamp,:]
groupinds = Utils.findgroups(task_vars)
groups = OrderedDict(k=>task_vars[groupinds.==k,:][1,:] for k in unique(groupinds))
cp_cond = OrderedDict()
pc_cond = OrderedDict()
for group in keys(groups)
    inds = groupinds.==group
    ca1 = embedding_rows[:ca1][inds]
    pfc = embedding_rows[:pfc][inds]
    cp_cond[groups[group]] = Threads.@spawn CausalityTools.predictive_asymmetry(ca1,
                                                         pfc, est,
                                                         params.horizon);
    pc_cond[groups[group]] = Threads.@spawn CausalityTools.predictive_asymmetry(pfc,
                                                         ca1, est,
                                                         params.horizon);

end

cp_manifold_pa = fetch(cp_manifold_pa)
pc_manifold_pa = fetch(pc_manifold_pa)
pc_cond = Dict(k=>fetch(v) for (k,v) in pc_cond)
cp_cond = Dict(k=>fetch(v) for (k,v) in cp_cond)

serialize(datadir("manifold","manifold_pa_cause$(ntopt_string((;params...)))"),
                  (;params..., est, cp_manifold_pa, pc_manifold_pa, cp_cond, pc_cond)
                 )


#  ================
# Task space projoections
#  ================
embeddingH
areas = (:ca1,:pfc)
embeddingH[:ca1] = umap((Rca1')[:,1:nsamp], 8)
embeddingH[:pfc] = umap((Rpfc')[:,1:nsamp], 8)

proj(x,y) = x * y' * (y * y') * y

proj(embedding[:ca1], task_vars)

