using DrWatson
include(scriptsdir("manifold","Umap_deserialize.jl"))

#  ================
# CAUSALITY AND MANFIOLD
#  ================

using Serialization
using Utils.namedtup: ntopt_string
using CausalityTools
using Entropies
using Munge.causal
using Plots

@info "prop missinG"
esttype = :binned
est, params = get_est_preset(esttype)
params = (;params..., horizon=1:250, thread=true)
params = (;params..., binning=5)

@assert !isempty(embedding)
CA1PFC = Munge.causal.get_paired_embeddings(embedding, :ca1, :pfc)
PFCCA1 = Munge.causal.get_paired_embeddings(embedding, :pfc, :ca1)

# Test works
#embeddingX,embeddingY = first(values(CA1PFC))
#CausalityTools.predictive_asymmetry(
#                                    eachrow(Dataset(embeddingX)), 
#                                    eachrow(Dataset(embeddingY)), 
#                                    est,
#                                    params[:horizon])

Threads.nthreads() = 16
Threads.nthreads()
using DiskBackedDicts, ThreadSafeDicts
global G_ca1pfc = G_pfcca1 = C_ca1pfc = C_pfcca1 = nothing
if params[:thread]
    Dtype = ThreadSafeDict
    G_ca1pfc, G_pfcca1 = ThreadSafeDict(), ThreadSafeDict()
    C_ca1pfc, C_pfcca1 = ThreadSafeDict(), ThreadSafeDict()
else
    Dtype = DiskBackedDict
    G_ca1pfc, G_pfcca1 = DiskBackedDict("G_ca1pfc.jld2"), DiskBackedDict("G_pfcca1.jld2")
    C_ca1pfc, C_pfcca1 = DiskBackedDict("C_ca1pfc.jld2"), DiskBackedDict("C_pfcca1.jld2")
end
GC.gc()

conditional_pred_asym(C_ca1pfc, CA1PFC, beh, [:cuemem,:correct]; 
                      groups=[[1,1],[1,0],[0,1],[0,0]], 
                      inds_of_t, params...)
conditional_pred_asym(C_pfcca1, PFCCA1, beh, [:cuemem,:correct]; 
                      groups=[[1,1],[1,0],[0,1],[0,0]], 
                      inds_of_t, params...)

predictive_asymmetry!(G_ca1pfc, CA1PFC; params...)
predictive_asymmetry!(G_pfcca1, PFCCA1; params...)

Idone=hcat([[!ismissing(vv) && istaskdone(vv) for vv in values(v)]   for v in [V for V in values(C_pfcca1)]]...)
Ifail=hcat([[!ismissing(vv) && istaskfailed(vv) for vv in values(v)] for v in [V for V in values(C_pfcca1)]]...)
Imissing=hcat([[!ismissing(vv) && istaskfailed(vv) for vv in values(v)] for v in [V for V in values(C_pfcca1)]]...)

condition_data_plot = plot(heatmap(Idone, title="done"), heatmap(Ifail,title="fail"), heatmap(Imissing,title="missing"))

# SAVING
paramstr = Utils.namedtup.tostring(params)
tagstr = if "tag" in propertynames(Main)
    "_$tag"
else
    "$animal$day.$(N)seg"
end
savefile = datadir("manifold","causal","pa_cause_$(paramstr)_$tagstr.jld2")
using JLD2
JLD2.@save(savefile, params, est, G_ca1pfc, G_pfcca1, C_ca1pfc, C_pfcca1,
          Idone, Ifail, Imissing)

using JLD2
D=JLD2.load(savefile)


"""

*=====*
|PLOTS|
*=====*

"""

#em = Table.to_dataframe(embedding, explode=false)
Plot.setfolder("manifold", "trust")
em = Munge.causal.make_embedding_df(embedding, inds_of_t, scores, beh)
histogram(em.score)
Plot.save(tagstr)
#opt = getoptions(em)

#Plot.setfolder("manifold","GEN_CAUSAL")
#plotcause(mean(values(filter(v-> v[1].n_neighbors==5, G_ca1pfc))), alpha=0.4, label="CA1->PFC");
#plotcause!(mean(values(filter(v->v[1].n_neighbors==5, G_pfcca1))), alpha=0.4, label="PFC->CA1")
#xlims!((0,250))
#Plot.save((;desc="mean",n_neighbors=5,min_dist=0.3))
#
#
#em=filter(v -> v[1].n_neighbors==5 && v[1].dim==3 && v[1].area==:ca1, embedding)
#sc1=scatter(eachcol(first(values(em)))..., markersize=1)
#em=filter(v -> v[1].n_neighbors==5 && v[1].dim==3 && v[1].area==:pfc, embedding)
#sc2=scatter(eachcol(first(values(em)))..., markersize=1)
#plot(sc1,sc2)
#Plot.save((;desc="manifolds",n_neighbors=5,min_dist=0.3))
#
#
#Plot.setfolder("manifold","GEN_CAUSAL")
#plotcause(mean(values(filter(v-> v[1].n_neighbors==n_neighbors, G_ca1pfc))), alpha=0.4, label="CA1->PFC");
#plotcause!(mean(values(filter(v->v[1].n_neighbors==n_neighbors, G_pfcca1))), alpha=0.4, label="PFC->CA1")
#xlims!((0,250))
#Plot.save((;desc="mean",n_neighbors=5,min_dist=0.3))
#
#em=filter(v -> v[1].n_neighbors==n_neighbors && v[1].dim==3 && v[1].area==:ca1, embedding)
#sc1=scatter(eachcol(first(values(em)))..., markersize=1)
#em=filter(v -> v[1].n_neighbors==n_neighbors && v[1].dim==3 && v[1].area==:pfc, embedding)
#sc2=scatter(eachcol(first(values(em)))..., markersize=1)
#plot(sc1,sc2)
#Plot.save((;desc="manifolds",n_neighbors=5,min_dist=0.3))

K = filter(k->k.min_dist ∈ min_dist && k.n_neighbors ∈ n_neighbors &&
           k.metric ∈ metric && k.dim == dim && k.feature == feature, keys(G_ca1pfc))

using Plot.cause
function getcausedistovertime(C::Dict)
    V = [V for V in values(C) if !ismissing(V)]
    [(i * 1/30,V[j][i]) for i in 1:250 for j in 1:length(V)]
end

function plotcausedistovertime(C::Dict;ch=:black,cmc=:black,labelmc="",histalpha=0.6,kws...)
    if K !== nothing
        C = Dict(k=>v for (k,v) in C if k ∈ K)
    end
    V = [V for V in values(C) if !ismissing(V)]
    histme = [(i * 1/30,V[j][i]) for i in 1:250 for j in 1:length(V)]
    histogram2d(histme;kws...,alpha=histalpha)
    plotmediancause!(V; timefact=1/30, c=cmc, markersize=1,label=labelmc,
                    fill=0, fillalpha=0.35, linewidth=0.2)
    hline!([0], c=ch, linewidth=3, linestyle=:dash, 
           xlabel="time (s)", ylabel="ℙ_asymmetry", label="")
end

function getmedian(set_cause)
    set_cause = collect(values(set_cause))
    inds = [isassigned(set_cause,p) && !ismissing(p) for p in eachindex(set_cause)]
    set_cause = skipmissing(set_cause[inds])
    set_cause = hcat(set_cause...)'
    vec(median(transform(set_cause),dims=2))
end
function getmean(set_cause)
    set_cause = collect(values(set_cause))
    inds = [isassigned(set_cause,p) && !ismissing(p) for p in eachindex(set_cause)]
    set_cause = skipmissing(set_cause[inds])
    set_cause = hcat(set_cause...)'
    vec(mean(transform(set_cause),dims=2))
end

# --------------------------------------------------

Plot.setfolder("manifold","GEN_CAUSAL")

caukws=(;bins=2 .* (30,45))
plot(plotcausedistovertime(G_ca1pfc; cmc=:red, labelmc="ca1 → pfc",caukws...),
     plotcausedistovertime(G_pfcca1; cmc=:blue, labelmc="pfc → ca1",caukws...),
     ylim=(-0.0010, 0.0010), 
     size=(1200,600)
   )

Plot.save("GEN_CAUSAL-$tagstr")

caukws=(;bins=2 .* (30,45))
plot(plotcausedistovertime(G_ca1pfc; cmc=:red, labelmc="ca1 → pfc",caukws...,histalpha=1),
     plotcausedistovertime(G_pfcca1; cmc=:blue, labelmc="pfc → ca1",caukws...,histalpha=1),
     ylim=(-0.0010, 0.0010), 
     size=(1200,600)
   )

Plot.save("GEN_CAUSAL-$tagstr-opaque")

# --------------------------------------------------

Plot.setfolder("manifold","COND_CAUSAL")

condkws = (ylim=(-0.003,0.003),size=(900,900))
P1= plot(
         plotcausedistovertime(C_ca1pfc[[1,1]];title="CA1→PFC, MEM correct",cmc=:red,caukws...,bins=(60,200)),
         plotcausedistovertime(C_pfcca1[[1,1]];title="PFC→CA1, MEM correct",cmc=:blue,caukws...,bins=(60,150)),
         plotcausedistovertime(C_ca1pfc[[1,0]];title="CA1→PFC, MEM error",cmc=:red,caukws...,bins=(60,450)),
         plotcausedistovertime(C_pfcca1[[1,0]];title="PFC→CA1, MEM error",cmc=:blue,caukws...,bins=(60,300));
    condkws...
)

Plot.save("MEM-$tagstr")

P1= plot(
         plotcausedistovertime(C_ca1pfc[[1,1]];title="CA1→PFC, MEM correct",cmc=:red,caukws...,bins=(60,200), histalpha=1),
         plotcausedistovertime(C_pfcca1[[1,1]];title="PFC→CA1, MEM correct",cmc=:blue,caukws...,bins=(60,150), histalpha=1),
         plotcausedistovertime(C_ca1pfc[[1,0]];title="CA1→PFC, MEM error",cmc=:red,caukws...,bins=(60,450), histalpha=1),
         plotcausedistovertime(C_pfcca1[[1,0]];title="PFC→CA1, MEM error",cmc=:blue,caukws...,bins=(60,300), histalpha=1);
    condkws...
)
Plot.save("MEM-$tagstr-opaque")



P2= plot(
     plotcausedistovertime(C_ca1pfc[[0,1]];title="CA1→PFC, CUE correct",cmc=:red,caukws...,bins=(60,100)),
     plotcausedistovertime(C_pfcca1[[0,1]];title="PFC→CA1, CUE correct",cmc=:blue,caukws...,bins=(60,100)),
     plotcausedistovertime(C_ca1pfc[[0,0]];title="CA1→PFC, CUE error",cmc=:red,caukws..., bins=(60,300)),
     plotcausedistovertime(C_pfcca1[[0,0]];title="PFC→CA1, CUE error",cmc=:blue,caukws..., bins=(60,300));
    condkws...
)

Plot.save("CUE-$tagstr")

P2= plot(
     plotcausedistovertime(C_ca1pfc[[0,1]];title="CA1→PFC, CUE correct",cmc=:red,caukws...,bins=(60,100), histalpha=1),
     plotcausedistovertime(C_pfcca1[[0,1]];title="PFC→CA1, CUE correct",cmc=:blue,caukws...,bins=(60,100), histalpha=1),
     plotcausedistovertime(C_ca1pfc[[0,0]];title="CA1→PFC, CUE error",cmc=:red,caukws..., bins=(60,300), histalpha=1),
     plotcausedistovertime(C_pfcca1[[0,0]];title="PFC→CA1, CUE error",cmc=:blue,caukws..., bins=(60,300), histalpha=1);
    condkws...
)

Plot.save("CUE-$tagstr-opaque")


lab = Dict([0,1]=>"CUE correct", [0,0]=>"CUE error", [1,1]=>"MEM correct", [1,0]=>"MEM error")
Cm_ca1pfc=Dict("CA1→PFC "*lab[k] => getmean(v) for (k,v) in C_ca1pfc)
Cm_pfcca1=Dict("CA1→PFC "*lab[k] => getmean(v) for (k,v) in C_pfcca1)

