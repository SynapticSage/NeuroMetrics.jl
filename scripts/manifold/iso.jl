using GLM
using  StatsBase
using  StatsPlots

#  ================
#  ISOLATED SPIKING & MANFIOLD
#  ================

# QUESTION : Are manifold jumps larger during iso spikes?

# iso spikes
beh.index = 1:size(beh,1)
Utils.filtreg.register(beh,spikes,on="time",transfer=["index"])
lfp = DI.load_lfp("RY16", 36, tet=5);
lfp.time = lfp.time .- DI.min_time_records[1]
lfp = Munge.lfp.annotate_cycles(lfp)
Munge.spiking.isolated(spikes, lfp, include_samples=false)
spikes.super_isolated = spikes.isolated  .&& spikes.phase .> 3.14

# QUESTION : Are the predicted  from isolated mean

I = combine(groupby(spikes,:index), :isolated=>mean, :super_isolated=>mean)[1:nsamp,:]
E = sqrt.(sum(diff(embedding[:ca1], dims=1).^2,dims=2))
E = [0; E]
I[!,:E] .= E


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
