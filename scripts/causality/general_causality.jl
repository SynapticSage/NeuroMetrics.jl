#  ================
# CAUSALITY AND MANFIOLD
#  ================

using DrWatson
using Serialization
using Utils.namedtup: ntopt_string
using CausalityTools
using Entropies

if !hasproperty(Main, :esttype)
@info "prop missinG"
    esttype = :binned
end

embedding_rows = Dict(k=>Dataset(v) for (k,v) in embedding);
@info esttype
if esttype == :binned
    params = (;bins = 9, horizon=1:3000)
    est = VisitationFrequency(RectangularBinning(params.bins))
    @info "Starting binned estimator, "
elseif esttype == :symbolic
    params = (m = 15, τ = 1)
    #params = (m = 100, τ = 4)
    est = SymbolicPermutation(;params...)
    params = (;params..., horizon=1:3000)
else
    @warn "not recog"
end


cp_manifold_pa = Threads.@spawn CausalityTools.predictive_asymmetry(embedding_rows[:ca1],
                                                     embedding_rows[:pfc], est,
                                                     params.horizon);
pc_manifold_pa = Threads.@spawn CausalityTools.predictive_asymmetry(embedding_rows[:pfc],
                                                     embedding_rows[:ca1], est,
                                                     params.horizon);

cp_manifold_pa = fetch(cp_manifold_pa)
pc_manifold_pa = fetch(pc_manifold_pa)
serialize(datadir("manifold","manifold_pa_cause$(ntopt_string((;params...)))"),
                  (;params..., est, cp_manifold_pa, pc_manifold_pa)
                 )

