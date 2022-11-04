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
if !hasproperty(Main, :esttype)
@info "prop missinG"
    esttype = :binned
end
est, params = get_est_preset(esttype)
params = (;params...,binning=4)

CA1PFC = Munge.causal.get_paired_embeddings(embedding, :ca1, :pfc)
PFCCA1 = Munge.causal.get_paired_embeddings(embedding, :pfc, :ca1)

# Test works
#embeddingX,embeddingY = first(values(CA1PFC))
#CausalityTools.predictive_asymmetry(
#                                    eachrow(Dataset(embeddingX)), 
#                                    eachrow(Dataset(embeddingY)), 
#                                    est,
#                                    params[:horizon])

Threads.nthreads() = 3
Threads.nthreads()
G_ca1pfc = global_predictive_asymmetry(CA1PFC; params...)
G_pfcca1 = global_predictive_asymmetry(PFCCA1; params...)

paramstr = Utils.namedtup.tostring(params)
tagstr   = tag == "" ? tag : "_$tag"
savefile = datadir("manifold","causal","pa_cause_$(paramstr)$tagstr.serial")
serialize(savefile, (;params..., est, G_ca1pfc, G_pfcca1))

#serialize(datadir("manifold","manifold_pa_cause$(ntopt_string((;params...)))"),
#                  (;params..., est, cp_manifold_pa, pc_manifold_pa))
#

Plot.setfolder("manifold","GEN_CAUSAL")
plotcause(mean(values(filter(v-> v[1].n_neighbors==5, G_ca1pfc))), alpha=0.4, label="CA1->PFC");
plotcause!(mean(values(filter(v->v[1].n_neighbors==5, G_pfcca1))), alpha=0.4, label="PFC->CA1")
xlims!((0,250))
Plot.save((;desc="mean",n_neighbors=5,min_dist=0.3))


em=filter(v -> v[1].n_neighbors==5 && v[1].dim==3 && v[1].area==:ca1, embedding)
sc1=scatter(eachcol(first(values(em)))..., markersize=1)
em=filter(v -> v[1].n_neighbors==5 && v[1].dim==3 && v[1].area==:pfc, embedding)
sc2=scatter(eachcol(first(values(em)))..., markersize=1)
plot(sc1,sc2)
Plot.save((;desc="manifolds",n_neighbors=5,min_dist=0.3))


Plot.setfolder("manifold","GEN_CAUSAL")
plotcause(mean(values(filter(v-> v[1].n_neighbors==n_neighbors, G_ca1pfc))), alpha=0.4, label="CA1->PFC");
plotcause!(mean(values(filter(v->v[1].n_neighbors==n_neighbors, G_pfcca1))), alpha=0.4, label="PFC->CA1")
xlims!((0,250))
Plot.save((;desc="mean",n_neighbors=5,min_dist=0.3))

em=filter(v -> v[1].n_neighbors==n_neighbors && v[1].dim==3 && v[1].area==:ca1, embedding)
sc1=scatter(eachcol(first(values(em)))..., markersize=1)
em=filter(v -> v[1].n_neighbors==n_neighbors && v[1].dim==3 && v[1].area==:pfc, embedding)
sc2=scatter(eachcol(first(values(em)))..., markersize=1)
plot(sc1,sc2)
Plot.save((;desc="manifolds",n_neighbors=5,min_dist=0.3))
