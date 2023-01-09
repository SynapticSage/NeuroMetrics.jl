using GoalFetchAnalysis
import Plot
using Utils.namedtup: ntopt_string
using Munge.causal, Munge.manifold

using DrWatson
using Infiltrator
using Serialization
using CausalityTools
using Entropies
using DataStructures: OrderedDict
using Plots
using DataFrames
using Statistics, NaNStatistics, HypothesisTests
using StatsPlots
using JLD2

## ----------
## CONSTANTS
## ----------
corerr,tsk,cortsk = Munge.behavior.cor, Munge.behavior.tsk, 
              Munge.behavior.cortsk

## ----------
## PARAMETERS
## ----------
animal, day, N, filt = "RY16", 36, 100, nothing
areas = (:ca1,:pfc)
distance = :many
feature_engineer = :many # many | nothing
esttype = :binned
est, params = get_est_preset(esttype)
params = (;params..., horizon=1:30, thread=false)
params = (;params..., binning=5, window=1.25)
manifold.load_manis_workspace(Main, animal, day; filt, 
                              areas, distance, feature_engineer, 
                              N)

## -----------------------------
## COMPUTE
## -----------------------------

# GET PAIRED EMBEDDINGS
# (each set of embeddings takes about 15 seconds)
@assert !isempty(embedding)
@time CA1PFC = Munge.causal.get_paired_embeddings(embedding, :ca1, :pfc)
@time PFCCA1 = Munge.causal.get_paired_embeddings(embedding, :pfc, :ca1)

savefile = get_alltimes_savefile(animal, day, N; params)
JLD2.save(savefile; CA1PFC, PFCCA1, animal, day, N, params)

#= CONDITIONALs ( I no longer perform conditional in this script)
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

conditional_pred_asym(C_ca1pfc, CA1PFC, beh, conditionals; 
                      groups=[[1,1],[1,0],[0,1],[0,0]], 
                      inds_of_t, params...)
conditional_pred_asym(C_pfcca1, PFCCA1, beh, conditionals; 
                      groups=[[1,1],[1,0],[0,1],[0,0]], 
                      inds_of_t, params...)

predictive_asymmetry!(G_ca1pfc, CA1PFC; params...)
predictive_asymmetry!(G_pfcca1, PFCCA1; params...)

Idone=hcat([[!ismissing(vv) && istaskdone(vv) for vv in values(v)]   for v in [V for V in values(C_pfcca1)]]...)
Ifail=hcat([[!ismissing(vv) && istaskfailed(vv) for vv in values(v)] for v in [V for V in values(C_pfcca1)]]...)
Imissing=hcat([[!ismissing(vv) && istaskfailed(vv) for vv in values(v)] for v in [V for V in values(C_pfcca1)]]...)

condition_data_plot = plot(heatmap(Idone, title="done"), heatmap(Ifail,title="fail"), heatmap(Imissing,title="missing"))

# SAVING
using JLD2
JLD2.@save(savefile, params, est, G_ca1pfc, G_pfcca1, C_ca1pfc, C_pfcca1,
          Idone, Ifail, Imissing)
=#


#em = Table.to_dataframe(embedding, explode=false)
#opt = getoptions(em)

