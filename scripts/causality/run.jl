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
function Base.fetch(X::AbstractDict)
    Dict(k=>(try;fetch(v);catch;missing;end) for (k,v) in X)
end

## ----------
## CONSTANTS
## ----------
corerr,tsk,cortsk = Munge.behavior.cor, Munge.behavior.tsk, 
              Munge.behavior.cortsk


## ----------
## PARAMETERS
## ----------
animal, day, N, filt = "RY22", 21, 100, nothing
areas = (:ca1,:pfc)
distance = :many
feature_engineer = :many # many | nothing
esttype = :binned
est, params = get_est_preset(esttype)
params = (;params..., horizon=1:30, thread=true)
params = (;params..., binning=5, window=1.25)
manifold.load_manis_workspace(Main, animal, day; filt, 
                              areas, distance, feature_engineer, 
                              N)
spikes, beh, ripples, cells  = Load.load(animal, day)
savefile = get_alltimes_savefile(animal, day, N; params)

## -----------------------------
## COMPUTE
## -----------------------------

# GET PAIRED EMBEDDINGS
# (each set of embeddings takes about 15 seconds)
@assert !isempty(embedding)
@time CA1PFC = Munge.causal.get_paired_embeddings(embedding, :ca1, :pfc)
@time PFCCA1 = Munge.causal.get_paired_embeddings(embedding, :pfc, :ca1)

JLD2.jldsave(savefile; CA1PFC, PFCCA1, animal, day, N, params)

# SETUP
begin
    # Threading
    Threads.nthreads() = 16
    Threads.nthreads()

    # Obtain store for results
    using DiskBackedDicts, ThreadSafeDicts
    if params[:thread]
        Dtype = ThreadSafeDict
    else
        Dtype = Dict
    end
    predasym = Dtype()
    predasym["alltimes"] = Dtype("ca1pfc"=>Dtype(), "pfcca1"=>Dtype())
    info = Dtype()
    GC.gc()
end

# COMPUTE ALL TIMES
begin
    predictiveasymmetry!(predasym["alltimes"]["ca1pfc"], CA1PFC; params...)
    predictiveasymmetry!(predasym["alltimes"]["pfcca1"], PFCCA1; params...)
    G_ca1pfc = fetch(G_ca1pfc)
    G_pfcca1 = fetch(G_pfcca1)
end

function runconditional(ca1pfc, pfcca1, conditionals)
    conditional_pred_asym(ca1pfc, CA1PFC, beh, conditionals; 
                          groups=[[1,1],[1,0],[0,1],[0,0]], 
                          inds_of_t, params...)
    conditional_pred_asym(pfcca1, PFCCA1, beh, conditionals; 
                          groups=[[1,1],[1,0],[0,1],[0,0]], 
                          inds_of_t, params...)
end

# Obtain conditional runs
begin
    PROPS = [[:cuemem, :correct],[:cuemem,:correct,:hajtrajnum]]
    for props in PROPS

        🔑 = "props=" * replace(string(props),":"=>""," "=>"")
        if 🔑 ∉ keys(predasym)
            predasym[🔑], info[🔑] = Dtype(), Dtype()
        end

        C_ca1pfc, C_pfcca1 = getindex.(predasym[🔑], ["ca1pfc","pfcca1"])
        runconditional(C_ca1pfc, C_pfcca1, props)

        # Describe
        Idone=hcat([[!ismissing(vv) && istaskdone(vv) for vv in values(v)]  
                    for v in [V for V in values(C_pfcca1)]]...)
        Ifail=hcat([[!ismissing(vv) && istaskfailed(vv) for vv in values(v)] 
                    for v in [V for V in values(C_pfcca1)]]...)
        Imissing=hcat([[!ismissing(vv) && istaskfailed(vv) for vv in values(v)] 
                    for v in [V for V in values(C_pfcca1)]]...)
        condition_data_plot = plot(heatmap(Idone, title="done"), 
                                   heatmap(Ifail,title="fail"), heatmap(Imissing,title="missing"))
        setindex!.([info[🔑]], [Idone, Ifail, Imissing], ["done","fail","missing"])
        Plot.save("$animal.$day.$N.$props")

        C_ca1pfc = fetch(C_ca1pfc)
        C_pfcca1 = fetch(C_pfcca1)

        # Checkpoint
        using JLD2
        S=JLD2.jldsave(savefile; CA1PFC, PFCCA1, animal, day, N, params, est, 
                       predasym, info)
    end
end

# SAVING

#em = Table.to_dataframe(embedding, explode=false)
#opt = getoptions(em)

