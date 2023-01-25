#  ==================================
# PLOT ALL TIMES 
#  ==================================
using DrWatson
using Infiltrator
using ThreadSafeDicts
using JLD2, Serialization, CausalityTools, Entropies
using DataFrames, DataFramesMeta
using Statistics, NaNStatistics, HypothesisTests
using DataStructures: OrderedDict
using Plots, StatsPlots, ColorSchemes

using GoalFetchAnalysis
using Utils.namedtup: ntopt_string
import Plot
using Munge.manifold, Munge.causal, Munge.triggering
using Utils.binning
using Munge.causal
import Munge
using Plot.cause
using Plot.cause: plotmeancause, plotmediancause, plotmedianplushist, 
                   plotcausediff, getdiff, getmean, getmedian, 
                   getcausedistovertime

# =============
# Control panel
# =============
opt = Dict(
    # What transformations on the data?
    :diff => false,
    :summaryget => getmean,
    :link => :y,  # :none, :x, :xy, :y
    # What plots?
    :ploton => true,
    :plot => Dict(
        :manfiold_checks => true,
        :gen_causal => true,
        :cond_causal => true,
        :cond_causal_earlintlate => false,
        :ha => true
    )
)
summaryget = opt[:summaryget]

arena_ca1pfc_color(i,n) = get(ColorSchemes.Blues, 0.15 + 0.85*(i/n))
arena_pfcca1_color(i,n) = get(ColorSchemes.Reds,  0.15 + 0.85*(i/n))
home_ca1pfc_color(i,n) = get(ColorSchemes.Greens, 0.15 + 0.85*(i/n))
home_pfcca1_color(i,n) = get(ColorSchemes.Oranges,  0.15 + 0.85*(i/n))

## ----------
## PARAMETERS
## ----------
animal,day, filt, N="RY22",21,nothing,100
animal,day, filt, N="RY16",36,nothing,100
spikes, beh, ripples, cells  = Load.load(animal, day)
areas = (:ca1,:pfc)
distance = :many
feature_engineer = :many # many | nothing
esttype = :binned
est, params = get_est_preset(esttype, horizon=1:60, thread=true, binning=7, window=1.25)
manifold.load_manis_workspace(Main, animal, day; filt, 
                              areas, distance, feature_engineer, 
                              N)
storage = load_alltimes_savefile(animal, day, N; params)
diffed, predasym = false, storage["predasym"]
if opt[:diff] && !diffed
    Plot.setprepend("diff_")
    predasym = getdiff(predasym)
end
opt[:ploton] ? Plot.on() : Plot.off()
Plot.setappend("$animal.$day.$N")

x_time = collect(1:last(params[:horizon])) .* 1/30

G_ca1pfc, G_pfcca1 = predasym["alltimes"]["ca1pfc"], predasym["alltimes"]["pfcca1"]
C_ca1pfc, C_pfcca1 = predasym["props=[cuemem,correct]"]["ca1pfc"], 
                     predasym["props=[cuemem,correct]"]["pfcca1"]
H_ca1pfc, H_pfcca1 = predasym["props=[cuemem,correct,hatraj]"]["ca1pfc"], 
                     predasym["props=[cuemem,correct,hatraj]"]["pfcca1"]
H_ca1pfc, H_pfcca1 = predasym["props=[cuemem,correct,hatraj]"]["ca1pfc"], 
                     predasym["props=[cuemem,correct,hatraj]"]["pfcca1"]
h_ca1pfc, h_pfcca1 =  predasym["props=[cuemem,correct,ha]"]["ca1pfc"], 
                      predasym["props=[cuemem,correct,ha]"]["pfcca1"]
Hm_ca1pfc, Hm_pfcca1 =  predasym["props=[cuemem,correct,hatraj,moving]"]["ca1pfc"], 
                      predasym["props=[cuemem,correct,hatraj,moving]"]["pfcca1"]

## ----------
## CONSTANTS
## ----------
corerr,tsk,lab = Munge.behavior.corerr, Munge.behavior.tsk, 
                 Munge.behavior.cortsk

## ----------
## PARAMETERS
## ----------
est, params = get_est_preset(:binned, binning=7, window=1.25, horizon=1:30, thread=false)
datasets = (("RY16", 36, 100), ("RY22", 21, 100))

(animal, day, N) = datasets[1]

## ----------
## LOAD
## ----------
spikes, beh, ripples, cells  = Load.load(animal, day)
storage = load_alltimes_savefile(animal, day, N; params)
tagstr = "$animal.$day.$N"


# ---------------
# Embedding trust
# ---------------
if opt[:plot][:manifold_check]
    # Can we trust the embedding?
    Plot.setfolder("manifold", "trust")
    em = Munge.causal.make_embedding_df(embedding, inds_of_t, scores, beh)
    histogram(em.score, bins=100)
    Plot.save(tagstr)

    # Embedding trust per area
    @df em histogram(:score, bins=100, grouping=:area)
    Plot.save(tagstr * "_split=area")

    @df em histogram(:score, bins=100, grouping=:feature)
    Plot.save(tagstr * "_split=feature")
end

# ------------------------------------------
# FUNCTIONS
# ------------------------------------------
K = filter(k->k.min_dist ∈ min_dist && k.n_neighbors ∈ n_neighbors &&
           k.metric ∈ metric && k.dim == dim && k.feature == feature, keys(G_ca1pfc))
K = nothing

# --------------------------------------------------
if opt[:plots][:gen_causal]
    Plot.setfolder("manifold","GEN_CAUSAL")

    caukws=(;bins=2 .* (30,100))
    plot(plotmedianplushist(G_ca1pfc; cmc=:red,  labelmc="ca1 → pfc",caukws...),
         plotmedianplushist(G_pfcca1; cmc=:blue, labelmc="pfc → ca1",caukws...);
         ylim=(-0.0010, 0.0010), 
         size=(1200,600),
         caukws...
    )

    Plot.save("GEN_CAUSAL-$tagstr")
    caukws=(;bins=2 .* (30,100))
    plot(plotmedianplushist(G_ca1pfc; cmc=:red, labelmc="ca1 → pfc",caukws...,histalpha=1),
         plotmedianplushist(G_pfcca1; cmc=:blue, labelmc="pfc → ca1",caukws...,histalpha=1);
         ylim=(-0.0010, 0.0010), 
         size=(1200,600),
         caukws
       )

    Plot.save("GEN_CAUSAL-$tagstr-opaque")
end

# ==================================================
# CONDITIONAL CAUSAL
# ==================================================
if opt[:plot][:cond_causal]

    Plot.setfolder("manifold","COND_CAUSAL")

    # --------------------------
    # A. Standard time course plots
    # --------------------------
    condkws = (ylim=(-0.003,0.003),size=(800,800))
    P1= plot(
             plotmedianplushist(C_ca1pfc[[1,1]];title="CA1→PFC, MEM correct",cmc=:red,caukws...,bins=(60,200)),
             plotmedianplushist(C_pfcca1[[1,1]];title="PFC→CA1, MEM correct",cmc=:blue,caukws...,bins=(60,150)),
             plotmedianplushist(C_ca1pfc[[1,0]];title="CA1→PFC, MEM error",  cmc=:red,caukws...,bins=(60,450)),
             plotmedianplushist(C_pfcca1[[1,0]];title="PFC→CA1, MEM error",  cmc=:blue,caukws...,bins=(60,300));
        condkws...
    )

    Plot.save("MEM-$tagstr")

    P1= plot(
             plotmedianplushist(C_ca1pfc[[1,1]];title="CA1→PFC, MEM correct",cmc=:red,caukws...,bins=(60,200), histalpha=1),
             plotmedianplushist(C_pfcca1[[1,1]];title="PFC→CA1, MEM correct",cmc=:blue,caukws...,bins=(60,150), histalpha=1),
             plotmedianplushist(C_ca1pfc[[1,0]];title="CA1→PFC, MEM error",cmc=:red,caukws...,bins=(60,450), histalpha=1),
             plotmedianplushist(C_pfcca1[[1,0]];title="PFC→CA1, MEM error",cmc=:blue,caukws...,bins=(60,300), histalpha=1);
        condkws...
    )

    Plot.save("MEM-$tagstr-opaque")



    P2= plot(
         plotmedianplushist(C_ca1pfc[[0,1]];title="CA1→PFC, CUE correct",cmc=:red,caukws...,bins=(60,100)),
         plotmedianplushist(C_pfcca1[[0,1]];title="PFC→CA1, CUE correct",cmc=:blue,caukws...,bins=(60,100)),
         plotmedianplushist(C_ca1pfc[[0,0]];title="CA1→PFC, CUE error",cmc=:red,caukws..., bins=(60,300)),
         plotmedianplushist(C_pfcca1[[0,0]];title="PFC→CA1, CUE error",cmc=:blue,caukws..., bins=(60,300));
        condkws...
    )

    Plot.save("CUE-$tagstr")

    P2= plot(
         plotmedianplushist(C_ca1pfc[[0,1]];title="CA1→PFC, CUE correct",cmc=:red,caukws...,bins=(60,100), histalpha=1),
         plotmedianplushist(C_pfcca1[[0,1]];title="PFC→CA1, CUE correct",cmc=:blue,caukws...,bins=(60,100), histalpha=1),
         plotmedianplushist(C_ca1pfc[[0,0]];title="CA1→PFC, CUE error",cmc=:red,caukws..., bins=(60,300), histalpha=1),
         plotmedianplushist(C_pfcca1[[0,0]];title="PFC→CA1, CUE error",cmc=:blue,caukws..., bins=(60,300), histalpha=1);
        condkws...
    )

    Plot.save("CUE-$tagstr-opaque")

    # Create shorcut structures
    C_ca1pfc=OrderedDict(k=>C_ca1pfc[k] for k in keys(lab))
    C_pfcca1=OrderedDict(k=>C_pfcca1[k] for k in keys(lab))
    Cd_ca1pfc=OrderedDict(k=>getdiff(C_ca1pfc[k]) for k in keys(lab))
    Cd_pfcca1=OrderedDict(k=>getdiff(C_pfcca1[k]) for k in keys(lab))
    Cm_ca1pfc=OrderedDict("CA1→PFC "*lab[k] => summaryget(v) for (k,v) in C_ca1pfc)
    Cm_pfcca1=OrderedDict("PFC→CA1 "*lab[k] => summaryget(v) for (k,v) in C_pfcca1)

    # --------------------------
    # B. Just the means
    # --------------------------
    # PLOT CA1 → PFC
    plot([plot(x_time,v; fill=0, label=k, c=:red) for (k,v) in Cm_ca1pfc]..., xlabel="time", ylabel="𝔸", alpha=0.5)
    # PLOT PFC → CA1
    plot([plot(x_time,v; fill=0, label=k, c=:skyblue) for (k,v) in Cm_pfcca1]..., xlabel="time",  ylabel="𝔸", alpha=0.5)


    # --------------------------
    # C1. Jacknived flows
    # --------------------------
    kws = (;link=opt[:link])
    Plot.setfolder("manifold","COND_CAUSAL", "jacknived flows")
    getjacknives(x) = leaveoneout(x; func=func_full)
    # How do those bins look on my image?
    plot([(plot(x_time, summaryget(v); fill=0, label=k, c=:skyblue); 
           plot!(x_time, getjacknives(v); c=:black, linestyle=:dash, alpha=0.25, label="");
           vline!(collect(Iterators.flatten(bins)), c=:black, linestyle=:dash, label=""))
           for (k,v) in C_ca1pfc]..., xlabel="time", ylabel="𝔸", alpha=0.5, kws...)
    Plot.save((;direction="ca1-pfc", params...))

    # PLOT PFC → CA1
    plot([(plot(x_time, summaryget(v); fill=0, label=k, c=:red); 
           plot!(x_time, getjacknives(v); c=:black, linestyle=:dash, alpha=0.25, label="");
           vline!(collect(Iterators.flatten(bins)), c=:black, linestyle=:dash, label=""))
          for (k,v) in C_pfcca1]..., xlabel="time", ylabel="𝔸", alpha=0.5, kws...,
         ylim=(-0.003, 0.002))
    Plot.save((;direction="pfc-ca1", params...))

    # BOTH 
    p=plot([(plot(x_time, summaryget(v1); fill=0, label="ca1→pfc",  legend_title="direction",
                  title=lab[k1], c=:red, alpha=0.5); 
           plot!(x_time, getjacknives(v1); c=:black, linestyle=:dash, alpha=0.25, label="");
           plot!(x_time, summaryget(v2); fill=0, label="pfc→ca1 ", legend_title="direction", title=lab[k2], c=:skyblue, alpha=0.5); 
           plot!(x_time, getjacknives(v2); c=:black, linestyle=:dash, alpha=0.25, label="");
           vline!(collect(Iterators.flatten(bins)), c=:black, linestyle=:dash, label=""))
          for ((k1,v1),(k2,v2)) in zip(C_ca1pfc, C_pfcca1)]..., 
           xlabel="time", ylabel="𝔸", alpha=0.5, kws..., ylim=(-0.0034, 0.0034), xlim=(0,3.5), size=(1000,400))

    Plot.save((;direction="both ca1-pfc pfc-ca1", params...))

    # --------------------
    # C2. Summary of jacknived
    # --------------------
    Plot.setfolder("manifold","COND_CAUSAL")
    Plot.setappend("$animal-$day")
    Plot.save("summary all on same axis")

    # ---------------------------------------
    # D. "Early, intermediate, late" analysis
    # ---------------------------------------
    if opt[:plot][:cond_causal_earlintlate]
        B = bin_the_curves(values(Cm_ca1pfc)...; bins, x_time)
        V = collect(skipmissing(values(first(values(C_ca1pfc)))))
        J = [leaveoneout(V; func=func_bin) for V in collect(values(C_ca1pfc))]
        J = [vcat(j...) for j in J]
        CIs = [(quantile(j[:,col],0.025), quantile(j[:,col], 0.975)) for j in J, col in 1:3]
        yerrors = [abs.(c.-b) for (c,b) in zip(CIs,B)]
        G = if trans == false
            groupedbar(replace.(collect(keys(Cm_ca1pfc))," "=>"\n"), B;
                       errorbar=yerrors, linewidth=2, label=("early","intermediate","late"),
                       legend=:none,
                       grid=false,minorgrid=false,
                       alpha=0.5, ylabel="Binned 𝔸")
        else
            manual_trans(x) = [x[i,j] for j in 1:size(x,2), i in 1:size(x,1)]
            groupedbar(["early","intermediate","late"], manual_trans(B);
                       legend=:none,
                       errorbar=yerrors, linewidth=2, 
                       label=(string.(collect(keys(Cm_pfcca1)))),
                       grid=false, minorgrid=false,
                       alpha=0.5, ylabel="Binned 𝔸")
        end
        Plot.save("ca1pfc - groupedbar")

        B = bin_the_curves(values(Cm_pfcca1)...; bins, x_time)
        V = collect(skipmissing(values(first(values(C_pfcca1)))))
        J = [leaveoneout(V; func=func_bin) for V in collect(values(C_pfcca1))]
        J = [vcat(j...) for j in J]
        CIs = [(quantile(j[:,col],0.025), quantile(j[:,col], 0.975)) for j in J, col in 1:3]
        yerrors = [abs.(c.-b) for (c,b) in zip(CIs,B)]
        G = if trans == false
            groupedbar(replace.(collect(keys(Cm_pfcca1))," "=>"\n"), B;
                       legend=:none,
                       grid=false,minorgrid=false,
                       errorbar=yerrors, linewidth=2, label=("early","intermediate","late"),
                       alpha=0.5, ylabel="Binned 𝔸")
        else
            manual_trans(x) = [x[i,j] for j in 1:size(x,2), i in 1:size(x,1)]
            groupedbar(["early","intermediate","late"], manual_trans(B);
                       legend=:none,
                       errorbar=yerrors, linewidth=2, 
                       grid=false,minorgrid=false,
                       label=(string.(collect(keys(Cm_pfcca1)))),
                       alpha=0.5, ylabel="Binned 𝔸")
        end
        Plot.save("pfcca1 - groupedbar")
    end



    # ------------
    # Summarizing the statsics of the flow
    # ------------
    v = first(values(C_pfcca1))
    @time Cb_pfcca1 = OrderedDict(k=>bin_the_curves(collect(values(v))...) 
                                  for (k,v) in C_pfcca1)
    @time Cb_ca1pfc = OrderedDict(k=>bin_the_curves(collect(values(v))...) 
                                  for (k,v) in C_ca1pfc)

    v = first(values(C_pfcca1))
    heatmap(v;clim=0.07 .* (-1,1), c=:vik, 
            title="visualizing the raw data that comes out")

    # Bootrap summaries


    # Hypothesis tests
    combos = Iterators.product( zip(["ca1pfc","pfcca1"], [Cb_ca1pfc, Cb_pfcca1]), 
                                Iterators.product([1,0],[1,0]))

    (((dir1, dict1), (task1, correct1)), ((dir2, dict2), (task2, correct2))) = first(Iterators.product(combos,combos))
    results = []
    @showprogress for (((dir1, dict1), (task1, correct1)),
                       ((dir2, dict2), (task2, correct2))) in Iterators.product(combos,combos)
        @info "loop" dir1 task1 correct1 dir2 task2 correct2
        try
        measures1 = dict1[[task1, correct1]]
        measures2 = dict2[[task2, correct2]]
        @softscope for (i,j) in Iterators.product(1:length(bins), 1:length(bins))
            @info "inner" i j
            unvarttest = HypothesisTests.UnequalVarianceTTest( measures1[:,i], measures2[:,j] )
            pval_ttest = pvalue(unvarttest)
            row = (;task1, task2, correct1, correct2, 
                   dir1,dir2, dir="$dir1 vs $dir2",
                   task="$(tsk[task1]) $(tsk[task2])",
                   correct="$(cor[correct1]) $(cor[correct2])",
                   unvarttest, pval_ttest
                  )
            row = OrderedDict(k=>[v] for (k,v) in zip(keys(row), values(row)))
            push!(results, DataFrame(row))
        end
        catch
            @infiltrate
        end
    end
    results = vcat(results...)
end

#                                       
# |   |,---.    |                  o    
# |---||---|    |--- ,---.,---.    .    
# |   ||   |    |    |    ,---|    |    
# `   '`   '    `---'`    `---^    |    
function get_ha_vals(set, f=nothing)
    sort([k for k in 
     collect(keys(H_ca1pfc))
     if k[1:length(set)] == set
     && (f===nothing || startswith(k[end],f))
    ])
end
function get_ha_vals(set)
    get_ha_vals(set, nothing)
end
#                              `---'    
if opt[:plot][:ha]

    kws = (;bins=(60,200), fillrange=0, fillalpha=0.025, linewidth=2)

    Plot.setfolder("causality","HATRAJ_CAUSAL")
    begin
        begin
            p1=plot()
            plot!(x_time, summaryget(H_ca1pfc[[0,1,"a1"]]);   title="CA1→PFC",  
                                          c=arena_ca1pfc_color(1,4),kws..., label="Cue-A1", linestyle=:dot)
            plot!(x_time, summaryget(H_ca1pfc[[0,1,"a2"]]);   
                                          c=arena_ca1pfc_color(2,4),kws..., label="Cue-A2", linestyle=:dot)
            plot!(x_time, summaryget(H_ca1pfc[[1,1,"a3"]]);  
                                          c=arena_ca1pfc_color(3,4),kws..., label="Mem-A1")
            plot!(x_time, summaryget(H_ca1pfc[[1,1,"a4"]]); 
                                          c=arena_ca1pfc_color(4,4),kws..., label="Mem-A2")
            #plot!(summaryfun(H_ca1pfc[[1,1,"a5"]]);   title="CA1→PFC, MEM correct",  
            #                              c=ca1pfc_color(3,4),kws..., label="3", linestyle=:dot)
            #plot!(summaryfun(H_ca1pfc[[1,1,"a6"]]);   title="CA1→PFC, MEM correct",  
            #                               c=ca1pfc_color(4,4),kws...,label="4", linestyle=:dot)
            hline!([0];c=:black,linestyle=:dot,label="")

            p2=plot()
            plot!(x_time, summaryget(H_pfcca1[[0,1,"a1"]]);   title="PFC→CA1",  
                                          c=arena_pfcca1_color(1,4),kws..., label="Cue-A1", linestyle=:dot)
            plot!(x_time, summaryget(H_pfcca1[[0,1,"a2"]]);                                          
                                          c=arena_pfcca1_color(2,4),kws..., label="Cue-A2", linestyle=:dot)
            plot!(x_time, summaryget(H_pfcca1[[1,1,"a3"]]);                                          
                                          c=arena_pfcca1_color(3,4),kws..., label="Mem-A1")
            plot!(x_time, summaryget(H_pfcca1[[1,1,"a4"]]);                                          
                                          c=arena_pfcca1_color(4,4),kws..., label="Mem-A2")
            #plot!(summaryfun(H_pfcca1[[1,1,"a5"]]);  
            #                              c=pfcca1_color(3,4),kws..., label="3", linestyle=:dot)
            #plot!(summaryfun(H_pfcca1[[1,1,"a6"]]); 
            #                               c=pfcca1_color(4,4),kws...,label="4", linestyle=:dot)
            hline!([0];c=:black,linestyle=:dot,label="")
            nothing
        end
        arena = plot(p1,p2, background_color=:seashell2)
        Plot.save("arena flow")

        begin
            p1=plot()
            plot!(x_time, summaryget(H_ca1pfc[[0,1,"h1"]]);   title="CA1→PFC",  
                                          c=home_ca1pfc_color(1,4),kws..., label="Cue-H1", linestyle=:dot)
            plot!(x_time, summaryget(H_ca1pfc[[0,1,"h2"]]);   
                                          c=home_ca1pfc_color(2,4),kws..., label="Cue-H2", linestyle=:dot)
            plot!(x_time, summaryget(H_ca1pfc[[1,1,"h3"]]);  
                                          c=home_ca1pfc_color(3,4),kws..., label="Mem-H1")
            plot!(x_time, summaryget(H_ca1pfc[[1,1,"h4"]]); 
                                          c=home_ca1pfc_color(4,4),kws..., label="Mem-H2")
            #plot!(summaryfun(H_ca1pfc[[1,1,"a5"]]);   title="CA1→PFC, MEM correct",  
            #                              c=ca1pfc_color(3,4),kws..., label="3", linestyle=:dot)
            #plot!(summaryfun(H_ca1pfc[[1,1,"a6"]]);   title="CA1→PFC, MEM correct",  
            #                               c=ca1pfc_color(4,4),kws...,label="4", linestyle=:dot)
            hline!([0];c=:black,linestyle=:dot,label="")

            p2=plot()
            plot!(x_time, summaryget(H_pfcca1[[0,1,"h1"]]);   title="PFC→CA1",  
                                          c=home_pfcca1_color(1,4),kws..., label="Cue-H1", linestyle=:dot)
            plot!(x_time, summaryget(H_pfcca1[[0,1,"h2"]]);   
                                          c=home_pfcca1_color(2,4),kws..., label="Cue-H2", linestyle=:dot)
            plot!(x_time, summaryget(H_pfcca1[[1,1,"h3"]]);   
                                          c=home_pfcca1_color(3,4),kws..., label="Mem-H1")
            plot!(x_time, summaryget(H_pfcca1[[1,1,"h4"]]);   
                                          c=home_pfcca1_color(4,4),kws..., label="Mem-H2")
            #plot!(summaryfun(H_pfcca1[[1,1,"a5"]]);  
            #                              c=pfcca1_color(3,4),kws..., label="3", linestyle=:dot)
            #plot!(summaryfun(H_pfcca1[[1,1,"a6"]]); 
            #                               c=pfcca1_color(4,4),kws...,label="4", linestyle=:dot)
            hline!([0];c=:black,linestyle=:dot,label="")
            nothing
        end
        home = plot(p1,p2, background_color=:seashell2)
        Plot.save("home flow")

        (;arena, home)
    end
    plot(home,arena, layout=Plots.grid(2,1), background_color=:seashell2) 


    Plot.setfolder("causality","HA_CAUSAL")
    begin
        p1=plot()
        plot!(x_time, summaryget(h_ca1pfc[[0,1,'H']]);   title="CA1→PFC",  
                                      c=home_ca1pfc_color(2,4),kws...,  label="Cue-H", linestyle=:dot)
        plot!(x_time, summaryget(h_ca1pfc[[0,1,'A']]);   
                                      c=arena_ca1pfc_color(2,4),kws..., label="Cue-A", linestyle=:dot)
        plot!(x_time, summaryget(h_ca1pfc[[1,1,'H']]);  
                                      c=home_ca1pfc_color(2,4),kws...,  label="Mem-H")
        plot!(x_time, summaryget(h_ca1pfc[[1,1,'A']]); 
                                      c=arena_ca1pfc_color(2,4),kws..., label="Mem-A")
        #plot!(summaryfun(h_ca1pfc[[1,1,"a5"]]);   title="CA1→PFC, MEM correct",  
        #                              c=ca1pfc_color(3,4),kws..., label="3", linestyle=:dot)
        #plot!(summaryfun(h_ca1pfc[[1,1,"a6"]]);   title="CA1→PFC, MEM correct",  
        #                               c=ca1pfc_color(4,4),kws...,label="4", linestyle=:dot)
        hline!([0];c=:black,linestyle=:dot,label="")

        p2=plot()
        plot!(x_time, summaryget(h_pfcca1[[0,1,'H']]);   title="PFC→CA1",  
                                      c=home_pfcca1_color(2,4),kws...,  label="Cue-H", linestyle=:dot)
        plot!(x_time, summaryget(h_pfcca1[[0,1,'A']]);   
                                      c=arena_pfcca1_color(2,4),kws..., label="Cue-A", linestyle=:dot)
        plot!(x_time, summaryget(h_pfcca1[[1,1,'H']]);   
                                      c=home_pfcca1_color(2,4),kws...,  label="Mem-H")
        plot!(x_time, summaryget(h_pfcca1[[1,1,'A']]);   
                                      c=arena_pfcca1_color(4,4),kws..., label="Mem-A")
        #plot!(summaryfun(h_pfcca1[[1,1,"a5"]]);  
        #                              c=pfcca1_color(3,4),kws..., label="3", linestyle=:dot)
    jk    #plot!(summaryfun(h_pfcca1[[1,1,"a6"]]); 
        #                               c=pfcca1_color(4,4),kws...,label="4", linestyle=:dot)
        hline!([0];c=:black,linestyle=:dot,label="")
        nothing
    end
    plot(p1, p2, background_color=:seashell2)
    Plot.save("home and arena flow summaries")


    # --------------------------
    # moving - home/arena - traj
    # --------------------------
    function plot_movhatraj(correct::Int=1)

        Plot.setfolder("causality","HATRAJ_moving_CAUSAL")
        background_color = if correct == -1
            Plot.setprepend("error_")
            :lightsalmon
        else
            Plot.setprepend("")
            :seashell2
        end

        c = correct == 1 ? "" : "*"
        arena_labels = ["Cue-A1", "Cue-A2", "Mem-A1", "Mem-A2"]
        home_labels = replace.(arena_labels, "A"=>"H")
        arena_hatraj = c .* "a" .* string.(collect(1:4))
        home_hatraj = replace.(arena_hatraj, "a"=>"h")
        println(arena_labels, home_labels, arena_hatraj, home_hatraj)
        
        # MOVING ARENA
        begin
            lab="Move"
            p1=plot()
            for (i,(l,h)) in enumerate(zip(arena_labels, arena_hatraj))
                cuemem = parse(Int,last(h)) > 2
                linestyle = cuemem ? :solid : :dot
                key = [cuemem, correct, h, true]
                if key in keys(Hm_ca1pfc)
                    y=summaryget(Hm_ca1pfc[key])
                    if y === missing; continue; end
                    plot!(x_time, y;   
                          title="$lab\nCA1→PFC", c=arena_ca1pfc_color(i,4),kws..., 
                          label=l, linestyle)
                end
            end
            hline!([0];c=:black,linestyle=:dot,label="")

            p2=plot()
            for (i,(l, h)) in enumerate(zip(arena_labels, arena_hatraj))
                cuemem = parse(Int,last(h)) > 2
                linestyle = cuemem ? :solid : :dot
                key = [cuemem,correct,h,true]
                if key in keys(Hm_pfcca1)
                    y=summaryget(Hm_pfcca1[key])
                    if y === missing; continue; end
                    plot!(x_time, y;   
                          title="$lab\nPFC→CA1", c=arena_pfcca1_color(i,4),kws..., 
                          label=l, linestyle)
                end
            end
            hline!([0];c=:black,linestyle=:dot,label="")
            nothing
        end
        movingarena = plot(p1,p2; background_color)
        Plot.save("moving arena flow")

        # MOVING HOME
        begin
            lab="Move"
            p1=plot()
            for (i,(l,h)) in enumerate(zip(home_labels, home_hatraj))
                cuemem = parse(Int,last(h)) > 2
                linestyle = cuemem ? :solid : :dot
                key = [cuemem, correct, h, true]
                y = summaryget(Hm_ca1pfc[key])
                if ismissing(y); continue; end
                if key in keys(Hm_ca1pfc)
                    plot!(x_time, y;   
                          title="$lab\nCA1→PFC", c=home_ca1pfc_color(i,4),kws..., 
                          label=l, linestyle)
                end
            end
            hline!([0];c=:black,linestyle=:dot)

            p2=plot()
            for (i,(l, h)) in enumerate(zip(home_labels, home_hatraj))
                cuemem = parse(Int,last(h)) > 2
                linestyle = cuemem ? :solid : :dot
                key = [cuemem,correct,h,true]
                y = summaryget(Hm_pfcca1[key])
                if ismissing(y); continue; end
                if key in keys(Hm_pfcca1)
                    plot!(x_time, y;   
                          title="$lab\nPFC→CA1", c=home_pfcca1_color(i,4),kws..., 
                          label=l, linestyle)
                end
            end
            hline!([0];c=:black,linestyle=:dot,label="")
            nothing
        end
        movinghome = plot(p1,p2; background_color)
        Plot.save("moving home flow")

        # ARENA STILL
        begin
            lab="Still"
            p1=plot()
            for (i,(l,h)) in enumerate(zip(arena_labels, arena_hatraj))
                cuemem = parse(Int,last(h)) > 2
                linestyle = cuemem ? :solid : :dot
                key = [cuemem,correct,h,false]
                y = summaryget(Hm_ca1pfc[key])
                if ismissing(y); continue; end
                if key in keys(Hm_ca1pfc)
                    plot!(x_time, y;   
                          title="$lab\nCA1→PFC", c=arena_ca1pfc_color(i,4),kws..., 
                          label=l, linestyle)
                end
            end
            hline!([0];c=:black,linestyle=:dot,label="")

            p2=plot()
            for (i,(l, h)) in enumerate(zip(arena_labels, arena_hatraj))
                cuemem = parse(Int,last(h)) > 2
                linestyle = cuemem ? :solid : :dot
                key = [cuemem,correct,h,false]
                y = summaryget(Hm_pfcca1[key])
                if ismissing(y); continue; end
                if key in keys(Hm_pfcca1)
                    plot!(x_time, y;   
                          title="$lab\nPFC→CA1", c=arena_pfcca1_color(i,4),kws..., 
                          label=l, linestyle)
                end
            end
            hline!([0];c=:black,linestyle=:dot,label="")
            nothing
        end
        stillarena = plot(p1,p2; background_color)
        Plot.save("still arena flow")

        # HOME STILL
        begin
            lab="Still"
            p1=plot()
            for (i,(l,h)) in enumerate(zip(home_labels, home_hatraj))
                cuemem = parse(Int,last(h)) > 2
                linestyle = cuemem ? :solid : :dot
                key = [cuemem,correct,h,false]
                y = summaryget(Hm_ca1pfc[key])
                if ismissing(y); continue; end
                if key in keys(Hm_ca1pfc)
                    plot!(x_time, y;   
                          title="$lab\nCA1→PFC", c=home_ca1pfc_color(i,4),kws..., 
                          label=l, linestyle)
                end
            end
            hline!([0];c=:black,linestyle=:dot,label="")

            p2=plot()
            for (i,(l, h)) in enumerate(zip(home_labels, home_hatraj))
                cuemem = parse(Int,last(h)) > 2
                linestyle = cuemem ? :solid : :dot
                key = [cuemem,correct,h,false]
                y = summaryget(Hm_pfcca1[key])
                if ismissing(y); continue; end
                if key in keys(Hm_pfcca1)
                    plot!(x_time, y;   
                          title="$lab\nPFC→CA1", c=home_pfcca1_color(i,4),kws..., 
                          label=l, linestyle)
                end
            end
            hline!([0];c=:black,linestyle=:dot,label="")
        end
        stillhome = plot(p1, p2; background_color)
        Plot.save("still home flow")

        alles = plot(movinghome, movingarena, stillhome, stillarena; size=(1200,800), layout=grid(2,2), background_color)
        Plot.save("move+still home+arena flow summaries")

        (;all=alles, movinghome, movingarena, stillhome, stillarena)

    end # plot_movhatraj (moving - home/arena - traj)

    C_movhatraj = plot_movhatraj(1)
    C_movhatraj.all

    I_movhatraj = plot_movhatraj(0)
    I_movhatraj.all
end
