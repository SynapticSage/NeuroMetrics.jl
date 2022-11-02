using DrWatson
#include(scriptsdir("manifold","Umap_deserialize.jl"))

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
    esttype = :symbolic
end
est, params = get_est_preset(esttype)

CA1PFC = Munge.causal.get_paired_embeddings(embedding, :ca1, :pfc)
PFCCA1 = Munge.causal.get_paired_embeddings(embedding, :pfc, :ca1)

# Test works
#embeddingX,embeddingY = first(values(CA1PFC))
#CausalityTools.predictive_asymmetry(
#                                    eachrow(Dataset(embeddingX)), 
#                                    eachrow(Dataset(embeddingY)), 
#                                    est,
#                                    params[:horizon])

Threads.nthreads() = 2
Threads.nthreads()
G_ca1pfc = global_predictive_asymmetry(CA1PFC, est; params...)
G_pfcca1 = global_predictive_asymmetry(PFCCA1, est; params...)

paramstr = Utils.namedtup.tostring(params)
savefile = datadir("manifol","causal","pa_cause_$(paramstr).serial")

#serialize(datadir("manifold","manifold_pa_cause$(ntopt_string((;params...)))"),
#                  (;params..., est, cp_manifold_pa, pc_manifold_pa))
#
