"""
Obtains the causality information in general for a whole ℳ(time)

where ℳ  is a manifold as a function of time

# Params
filt selects a condition

# Assumptions
you already loaded some data we can work on
(see causality/run.jl)
"""

using DrWatson
using Serialization
using Utils.namedtup: ntopt_string
using CausalityTools
using Entropies
using ProgressMeter

if !hasproperty(Main, :esttype)
@info "prop missinG"
    esttype = :binned
end

embedding_rows = Dict(k=>size(v,1)>size(v,2) ? Dataset(v) : Dataset(v')
                      for (k,v) in embedding);

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

if filt !== nothing
    params = (;params..., filt)
end

spikes, beh, ripples, cells  = Load.load(animal, day)
task_vars = replace([beh.cuemem beh.correct],NaN=>-1)[1:nsamp, :]
groupinds = Utils.findgroups(task_vars)
groups    = OrderedDict(k=>task_vars[groupinds.==k,:][1,:]
                      for k in unique(groupinds))
groups    = [(cuemem=1,correct=1),(cuemem=0,correct=1),(cuemem=1,correct=0),(cuemem=0,correct=0)]

cp_cond = OrderedDict()
pc_cond = OrderedDict()

@showprogress for group in keys(groups)
    inds = groupinds .== group
    ca1 = embedding_rows[:ca1][inds]
    pfc = embedding_rows[:pfc][inds]
    cp_cond[groups[group]] = Threads.@spawn CausalityTools.predictive_asymmetry(ca1,
                                                         pfc, est,
                                                         params.horizon);
    pc_cond[groups[group]] = Threads.@spawn CausalityTools.predictive_asymmetry(pfc,
                                                         ca1, est,
                                                         params.horizon);
end

pc_cond = Dict(k=>(try fetch(v); catch; end) for (k,v) in pc_cond)
cp_cond = Dict(k=>(try fetch(v); catch; end) for (k,v) in cp_cond)

serialize(datadir("manifold",
                  "manifold_task_pa_cause$(ntopt_string((;params...)))"),
                  (;params..., est, cp_cond, pc_cond)
                 )



plotcause(cp_cond[[1,  1]]; xlim=(0, 15), marker='.', alpha=0.2, label="c->p", rect=false)
plotcause!(pc_cond[[1, 1]], xlim=(0, 15), marker='.', alpha=0.2, label="p->c", rect=false)

plotcause(cp_cond[[0,  1]], xlim=(0, 15), marker='.', alpha=0.2, label="c->p")
plotcause!(pc_cond[[0, 1]], xlim=(0, 15), marker='.', alpha=0.2, label="p->c")

plotcause(cp_cond[[1,  0]], xlim=(0, 15), marker='.', alpha=0.2, label="c->p")
plotcause!(pc_cond[[1, 0]], xlim=(0, 15), marker='.', alpha=0.2, label="p->c")

plotcause(cp_cond[[0,  0]], xlim=(0, 15), marker='.', alpha=0.2, label="c->p")
plotcause!(pc_cond[[0, 0]], xlim=(0, 15), marker='.', alpha=0.2, label="p->c")
