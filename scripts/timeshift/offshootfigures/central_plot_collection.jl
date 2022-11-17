using GoalFetchAnalysis
using Serialization
using Plots
using Timeshift.shiftmetrics
using Field.metrics
using Timeshift
using DimensionalData
using Statistics, NaNStatistics, Bootstrap
import Plot
using StatsBase
using ProgressMeter
using DataFramesMeta
using SoftGlobalScope
theme(:bright)

shifts=-2.0f0:0.05f0:2.0f0
savefile = "/home/ryoung/Projects/goal-code/data/timeshift/fixed_shifts_$shifts.serial"
F = load_fields()
#F,I,shifts = deserialize(savefile);


datasets = 
(
 ("RY16", 36, "CA1"),
 ("RY22", 21, "CA1"),
 ("RY16", 36, "PFC"),
 ("RY22", 21, "PFC"),
 ("RY16", 36, nothing),
 ("RY22", 21, nothing),
 ("super", 0, "CA1"),
 ("super", 0, "PFC"),
 ("super", 0, nothing),
)

# If this doesn't work subfunctions may not be scoped right and may need SoftGlobalScope
animal, day, brain_area = datasets[end]
#include(expanduser("~/tmp2.jl"))

@softscope for (i,(animal,day,brain_area)) in collect(enumerate(datasets))

    @info "loop" i animal day brain_area
    if brain_area === nothing
        @info "skip"
        continue
    end

    # = PREPLOT ACTIONS =#
    cells = Load.load_cells(animal, day)

    keyz = Dict(key.datacut=>key for key in filter(k->k.animal==animal && 
                                                   k.day == day, keys(F)))

    @eval Plot parent_folder = ["timeshift", "central_plot_collection.jl"]
    Plot.setappend((;animal,day,brain_area))

    function prep(f)    
        push_dims!(f)
        push_celltable!(f, cells)
        pop_metric!(f, :unit)
        push_metric!(f, metrics.bitsperspike)
        push_metric!(f, metrics.totalcount)
        push_metric!(f, metrics.maxrate)
        push_shiftmetric!(f, best_tau!; metric=:bitsperspike)
        push_metric!(f, metrics.bitsperspikeold)
        push_shiftmetric!(f, best_tau!; metric=:bitsperspikeold)
        f = f[vec(all(f[:totalcount] .> 50, dims=2) .&&
                  any(f[:bitsperspike] .> 0.5, dims=2)) ,:]
        if brain_area !== nothing
            f = f[vec(all(f[:area] .== brain_area, dims=2)), :]
        end
        f
    end
      
    function doheat(f,key)
        inds = sortperm(f[:bestshift_bitsperspike][:,1])
        bps  = f[:bitsperspike][inds, :]
        XX = hcat([Utils.norm_extrema(b) for b in eachrow(bps)]...)'
        heatmap(shifts, collect(1:size(XX,1)), XX, title="$(key.datacut)", size=2 .* (600,800))
    end
    # = --------------------------------------------------------------- =#


    """
    This tracks changes in the tau best shift
    """

    bins = 25
    Plot.setfolder("histograms, best shift")
    for k in (:all, :nontask, :cue, :mem_error, :cue_error, :task,:memory, :mem_correct, :cue_correct)
        key = keyz[k]
        k  = string(k)
        f   = prep(matrixform(F[key]))
        @info "$key" length(unique(f[:unit]))

        vals = f[shift=At(0)][:bestshift_bitsperspike]
        edges = LinRange(minimum(f[:shift]), maximum(f[:shift]), bins+1)
        xlim = (edges[begin], edges[end])

        #bl = bayesian_blocks(vals, prior=AIC(), resolution=1000)
        #support, density = to_pdf(bl)
        #plot(support, density)

        histogram(vals; bins, xlim, title=k, label="")
        #histogram(bl, bins=bl.edges, xlim, title=k, label="")

        Plot.save(k)
    end

bins=61
k=:all
Plot.setfolder("jerry-rigged superanimal histogram ")
H=[]
@showprogress for k in (:all, :nontask, :cue, :mem_error, :cue_error, :task,:memory, :mem_correct, :cue_correct)
    key1, key2 = keyz1[k], keyz2[k]
    k  = string(k)
    f1,f2   = prep_ry16(matrixform(F[key1])), prep_ry22(matrixform(F[key2]))
    @info "$key" #"$(length(unique(f[:unit])))"

    vals1,vals2 = f1[shift=At(0)][:bestshift_bitsperspike], f2[shift=At(0)][:bestshift_bitsperspike]
    edges = LinRange(minimum(f[:shift]), maximum(f[:shift]), bins+1)
    xlim = (edges[begin], edges[end])

    #bl = bayesian_blocks(vals, prior=AIC(), resolution=1000)
    #support, density = to_pdf(bl)
    #plot(support, density)

   h= histogram(vcat(vals1, vals2); bins, xlim, title=k, label="")
   push!(H,h)
    vline!([0],c=:black,linestyle=:dash)
    Plot.save(k)
end

bins=61
k=:all
Plot.setfolder("jerry-rigged superanimal histogram  comparison")
Plot.setappend("animal=jerryrig")
H=[]
@showprogress for (k1,k2) in ((:cue,), (:memory,))

    key11, key12 = keyz1[k1], keyz2[k1]
    key21, key22 = keyz1[k2], keyz2[k2]

    k1, k2  = string(k1),string(k2)
    f11 = prep_ry16(matrixform(F[key11]))
    f12 = prep_ry22(matrixform(F[key12]))
    f21 = prep_ry16(matrixform(F[key21]))
    f22 = prep_ry22(matrixform(F[key22]))

    #@info "$key" #"$(length(unique(f[:unit])))"

    vals11,vals12,vals21,vals22 = f11[shift=At(0)][:bestshift_bitsperspike], 
                    f12[shift=At(0)][:bestshift_bitsperspike],
                    f21[shift=At(0)][:bestshift_bitsperspike],
                    f22[shift=At(0)][:bestshift_bitsperspike]

    edges = LinRange(minimum(f[:shift]), maximum(f[:shift]), bins+1)
    xlim = (edges[begin], edges[end])

   #bl = bayesian_blocks(vals, prior=AIC(), resolution=1000)
   #support, density = to_pdf(bl)
   #plot(support, density)

   V1=vcat(vals11, vals12)
   V2=vcat(vals21, vals22)
  h1 = histogram(V1;  bins=9, xlim, title=k1 * " " * k2, label=k1, alpha=0.5)
  h2 = histogram!(V2; bins=9, xlim, label=k2, alpha=0.5)
  Plot.save((;k1,k2))
  using HypothesisTests
  e1, e2 = ecdf(V1), ecdf(V2)
  K = KruskalWallisTest(e1.sorted_values, e2.sorted_values)
  plot(e1,label=k1)
  plot!(e2, title="pval=$(pvalue(K))", xlabel="best shift on bits per spike", label=k2)
  Plot.save("bits per spike BEST SHIFT IS NOT sig different")


  push!(H,h)
   vline!([0],c=:black,linestyle=:dash)
   Plot.save(k)
nd




   bins = 25
   Plot.setfolder("tracking change in τ best_shift")
   for (k1, k2) in zip((:nontask, :cue, :mem_error, :cue_error),
                       (:task,:memory, :mem_correct, :cue_correct))
   @info "keys" k1 k2
       key1,key2 = keyz[k1],keyz[k2]
       k1s, k2s = replace.(string.([k1,k2]),"_" => "\\;")
       f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]));
       u = intersect(f1.dims[1], f2.dims[1])

       vals = f1[unit=At(u)][:bestshift_bitsperspike],
              f2[unit=At(u)][:bestshift_bitsperspike]
              
        edges = LinRange(minimum(f1[:shift]), maximum(f1[:shift]), bins+1)
        xlim = (edges[begin], edges[end])
        
        shiftdist = (vals[2] - vals[1])[:,1]

        # HISTOGRAM ("PDF")
        plot_histdist = 
            histogram(shiftdist; bins=edges,  xlim, normalize=:probability,
                      alpha=0.5, label="",
                      xlabel="shift\$_{$k2s}\$-shift\$_{$k1s}\$",
                      ylabel="pdf", title="$k2 - $k1")
        vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
            
        # HISTOGRAM ("CDF")
        dist = fit(Histogram, shiftdist, edges)
        dist = StatsBase.normalize(dist, mode=:probability)
        plot_histdist_cdf = 
            bar(dist.edges, cumsum(dist.weights); xlim, alpha=0.5, label="",
                xlabel="shift", ylabel="cdf")
            vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
            
        future = sum(shiftdist .> 0 .&& abs.(shiftdist) .< 1.5)
        past   = sum(shiftdist .< 0 .&& abs.(shiftdist) .< 1.5)
        b_fp =   bar(["past","future"],[past,future]./(past+future),
                     legend=:none)
        locality_thresh = 0.1;
        locl     = sum(abs.(shiftdist) .<= locality_thresh)
        nonlocal = sum(abs.(shiftdist) .> locality_thresh)
        b_lc = bar(["proximal","distal"], [locl, nonlocal])

        b_ratios = bar(["future/past", "shift/remain"],
                       [future/past, nonlocal/locl]);
                        
        plot(plot_histdist, plot_histdist_cdf,
             b_fp,b_lc,
              layout=[grid(2,1) grid(2,1)])
              
        Plot.save((;k1,k2))
   end

   statfunc = mean
   Plot.setfolder("comparison of bps via statfunc=$statfunc")
   P=[]
   for (k1, k2) in zip((:nontask, :cue, :mem_error, :cue_error, :error),
                       (:task,:memory, :mem_correct, :cue_correct, :correct))
       @info "keys" k1 k2
       key1,key2 = keyz[k1],keyz[k2]
       f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]))
       u = intersect(f1.dims[1], f2.dims[1])
       f1,f2 = f1[unit=At(u)], f2[unit=At(u)]

       sort_bs = sortperm(f1[:bestshift_bitsperspike][:,1])
       XX1 =bps  = f1[:bitsperspike][sort_bs, :]
       V1 = vec(statfunc(XX1,dims=1))
       bs1 = bootstrap(statfunc, collect(eachrow(XX1)), BalancedSampling(1000))
       straps1 = hcat(straps(bs1)...)

       sort_bs = sortperm(f2[:bestshift_bitsperspike][:,1])
       XX2 = bps  = f2[:bitsperspike][sort_bs, :]
       V2 = vec(statfunc(XX2,dims=1))
       bs2 = bootstrap(statfunc, collect(eachrow(XX2)), BasicSampling(1000))
       straps2 = hcat(straps(bs1)...)

       # I can do this because I intersected neurons in both
       bs21 = bootstrap(statfunc, collect(eachrow(XX2-XX1)), BasicSampling(1000))
       ci21 = confint(bs21, BasicConfInt(0.95))
       straps21 = [x-y for (x,y) in Iterators.product(eachrow(straps1),
                                                      eachrow(straps2))]
       shifts = collect(f1.dims[2])
       # Manual jackknife
       M = Matrix(XX2-XX1)
       jacknives = [mean(M[setdiff(1:size(M,1),[i]),:],dims=1)
                    for i in 1:size(M,1)]
       plot(hcat([[y for y in x] for x in ci21]...)')
       
       p=plot(shifts, V2-V1, title="$k2 - $k1", size=0.75 .* (800,1200),
              fillrange=0, xlabel="shift", ylabel="percent bits per spike")
       vline!([0],c=:black,linestyle=:dash, linewidth=2)
       plot!(shifts,hcat([vec(m) for m in
                          jacknives]...),c=:black,legend=:none,
             linestyle=:dot, alpha=0.4,title="Jacknifed
             $(string(Symbol(statfunc))) difference in\n bits per spike,
             $k2-$k1");
       hline!([0],c=:red,linestyle=:dash)
       
       if (k1,k2) == (:cue, :memory)
           shifts_local = Utils.in_range(shifts, [-.21, .61])
           shifts_distal = (!).(shifts_local)
           get_local(val)  = nanmean(val[shifts_local])
           get_distal(val) = nanmean(val[shifts_distal])
           L, D = get_local((V2-V1)),
                  get_distal((V2-V1))
          stde(x) = nanstd(x,dims=1)/sqrt(size(x,1))
          J = vcat(jacknives...)
          Lse, Dse = mean(get_local(stde(J))), mean(get_distal(stde(J)))
          b=bar(["local", "nonlocal"],   [L,D], yerror=[Lse, Dse], label="")
          ##scatter!(["local", "nonlocal"], [L,D], err=[Lse, Dse],
                   #markerstrokecolor = :transparent,
                   #linestyle=:dash, c=:black, markersize=10, label="", lw=4)
          plot!(["local", "nonlocal"], [L,D], err=[Lse, Dse],
                   l=nothing, m=(4,stroke(:red, 1.5)), label="")

          
          plot(p, b, size=(1000,1000))
          Plot.save((;desc="barsum", k1, k2))
       end
       push!(P,p)
       p
       Plot.save((;xlim="full", k1, k2))
   end

   plot(P..., size=(1600,600), fillalpha=0.5)
   Plot.save("together")


   statfunc = mean
   Plot.setfolder("comparison of norm_bps via statfunc=$statfunc")
   P=[]
   for (k1, k2) in zip((:nontask, :cue, :mem_error, :cue_error, :error),
                       (:task,:memory, :mem_correct, :cue_correct, :correct))
       @info "keys" k1 k2
       key1,key2 = keyz[k1],keyz[k2]
       f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]))
       u = intersect(f1.dims[1], f2.dims[1])
       f1,f2 = f1[unit=At(u)], f2[unit=At(u)]

       sort_bs = sortperm(f1[:bestshift_bitsperspike][:,1])
       XX1 =bps  = vcat([Utils.norm_extrema(f)[Utils.na, :]
                    for f in 
                    eachrow(f1[:bitsperspike][sort_bs, :])]...)

       V1 = vec(statfunc(XX1,dims=1))
       bs1 = bootstrap(statfunc, collect(eachrow(XX1)), BalancedSampling(1000))
       straps1 = hcat(straps(bs1)...)

       sort_bs = sortperm(f2[:bestshift_bitsperspike][:,1])
       XX2 =bps  = vcat([Utils.norm_extrema(f)[Utils.na, :]
                    for f in 
                    eachrow(f2[:bitsperspike][sort_bs, :])]...)
       V2 = vec(statfunc(XX2,dims=1))
       bs2 = bootstrap(statfunc, collect(eachrow(XX2)), BasicSampling(1000))
       straps2 = hcat(straps(bs1)...)

       # I can do this because I intersected neurons in both
       bs21 = bootstrap(statfunc, collect(eachrow(XX2-XX1)), BasicSampling(1000))
       ci21 = confint(bs21, BasicConfInt(0.95))
       straps21 = [x-y for (x,y) in Iterators.product(eachrow(straps1),
                                                      eachrow(straps2))]
       shifts = collect(f1.dims[2])
       # Manual jackknife
       M = Matrix(XX2-XX1)
       jacknives = [mean(M[setdiff(1:size(M,1),[i]),:],dims=1)
                    for i in 1:size(M,1)]
       plot(hcat([[y for y in x] for x in ci21]...)')
       
       p=plot(shifts, V2-V1, title="$k2 - $k1", size=0.75 .* (800,1200), fillrange=0,
             xlabel="shift", ylabel="Δ percent bits per spike")
       vline!([0],c=:black,linestyle=:dash, linewidth=2)
       plot!(shifts,hcat([vec(m) for m in
                          jacknives]...),c=:black,legend=:none,
             linestyle=:dot, alpha=0.4,title="Jacknifed
             $(string(Symbol(statfunc))) difference in\n bits per spike,
             $k2-$k1");
       hline!([0],c=:red,linestyle=:dash)
       
       if (k1,k2) == (:cue, :memory)
           shifts_local = Utils.in_range(shifts, [-.21, .61])
           shifts_distal = (!).(shifts_local)
           get_local(val)  = nanmean(val[shifts_local])
           get_distal(val) = nanmean(val[shifts_distal])
           L, D = get_local((V2-V1)),
                  get_distal((V2-V1))
          stde(x) = nanstd(x,dims=1)/sqrt(size(x,1))
          J = vcat(jacknives...)
          Lse, Dse = mean(get_local(stde(J))), mean(get_distal(stde(J)))
          b=bar(["local", "nonlocal"],   [L,D], yerror=[Lse, Dse], label="")
          ##scatter!(["local", "nonlocal"], [L,D], err=[Lse, Dse],
                   #markerstrokecolor = :transparent,
                   #linestyle=:dash, c=:black, markersize=10, label="", lw=4)
          plot!(["local", "nonlocal"], [L,D], err=[Lse, Dse],
                   l=nothing, m=(4,stroke(:red, 1.5)), label="")

          
          plot(p, b, size=(1000,1000))
          Plot.save((;desc="barsum", k1, k2))
       end
       push!(P,p)
       p
       Plot.save((;xlim="full", k1, k2))
   end

   plot(P..., size=(1600,600), fillalpha=0.5)
   Plot.save("together")

   """
   -------------------------------------------------------------------
   Obtain delta heatmaps
   """
   Plot.setfolder("delta in bps heatmap - difference index")
   for (k1, k2) in zip((:nontask, :cue, :correct, :mem_correct),
                       (:task,:memory, :error, :mem_error),)
       key1,key2 = keyz[k1],keyz[k2]
       f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]))
       D = intersect(f1.dims[1],f2.dims[1])
       R, L = [], []
       for d in D
           val = (f1[unit=At(d)][:bitsperspike]-f2[unit=At(d)][:bitsperspike])./
                     (f1[unit=At(d)][:bitsperspike]+f2[unit=At(d)][:bitsperspike])
           push!(R,val)
           push!(L, argmax(val))
       end
        ind = sortperm(L)
        R = hcat(R...)'
        lay = @layout [a; b{0.2h}]
        plot(heatmap(R[ind,:], clim=(-0.5,0.5)), plot(mean(R,dims=1)'),
             layout=lay)
        Plot.save((;k1,k2))
   end

   """
   -------------------------------------------------------------------
   """

   Plot.setfolder("delta in bps heatmap - difference index - norm01")
   for (k1, k2) in zip((:task,:memory, :error, :mem_error), (:nontask, :cue,
                                                             :correct,
                                                             :mem_correct),)
       key1,key2 = keyz[k1],keyz[k2]
       f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]))
       D = intersect(f1.dims[1],f2.dims[1])
       R, L = [], []
       for d in D
           val = (f2[unit=At(d)][:bitsperspike]-f1[unit=At(d)][:bitsperspike])./
               (f2[unit=At(d)][:bitsperspike]+f1[unit=At(d)][:bitsperspike])
           val = Utils.norm_extrema(val,[0,1])
           push!(R,val)
           push!(L, argmax(val))
       end
        ind = sortperm(L)
        R = hcat(R...)'
        lay = @layout [a; b{0.2h}]
        plot(heatmap(R[ind,:], clim=(0,1),c=:vik), plot(mean(R,dims=1)'),
             layout=lay)
        Plot.save((;k1,k2))
   end

   Plot.setfolder("delta in bps heatmap -- raw")
   for (k1, k2) in zip((:task,:memory, :error, :mem_error), (:nontask, :cue,
                                                             :correct,
                                                             :mem_correct),)
       key1,key2 = keyz[k1],keyz[k2]
       f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]))
       D = intersect(f1.dims[1],f2.dims[1])
       R, L = [], []
       for d in D
           val = f2[unit=At(d)][:bitsperspike]-f1[unit=At(d)][:bitsperspike]
           val = Utils.norm_extrema(val,[0,1])
           push!(R,val)
           push!(L, argmax(val))
       end
        ind = sortperm(L)
        R = hcat(R...)'
        lay = @layout [a; b{0.2h}]
        plot(heatmap(R[ind,:], clim=(0,1),c=:vik), plot(mean(R,dims=1)'),
             layout=lay)
        Plot.save((;k1,k2))
   end

   Plot.setfolder("delta in bps heatmap - difference index - norm01")
   for (k1, k2) in zip((:task,:memory, :error, :mem_error), (:nontask, :cue,
                                                             :correct,
                                                             :mem_correct),)
       key1,key2 = keyz[k1],keyz[k2]
       f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]))
       D = intersect(f1.dims[1],f2.dims[1])
       R, L = [], []
       for d in D
           val = (f2[unit=At(d)][:bitsperspike]-f1[unit=At(d)][:bitsperspike])./
               (f2[unit=At(d)][:bitsperspike]+f1[unit=At(d)][:bitsperspike])
           val = Utils.norm_extrema(val,[0,1])
           push!(R,val)
           push!(L, argmax(val))
       end
        ind = sortperm(L)
        R = hcat(R...)'
        lay = @layout [a; b{0.2h}]
        plot(heatmap(R[ind,:], clim=(0,1),c=:vik), plot(mean(R,dims=1)'),
             layout=lay)
        Plot.save((;k1,k2))
   end

   Plot.setfolder("delta in bps heatmap - norm01 => diff => norm01")
   for (k1, k2) in zip((:task,:memory, :error, :mem_error), (:nontask, :cue,
                                                             :correct,
                                                             :mem_correct),)
       key1,key2 = keyz[k1],keyz[k2]
       f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]))
       D = intersect(f1.dims[1],f2.dims[1])
       R, L = [], []
       for d in D
           b1, b2 = f1[unit=At(d)][:bitsperspike], f2[unit=At(d)][:bitsperspike]
           b1, b2 = Utils.norm_extrema(b1,[0,1]), Utils.norm_extrema(b2,[0,1])
           val = Utils.norm_extrema(b2-b1,[0,1])
           push!(R,val)
           push!(L, argmax(val))
       end
        ind = sortperm(L)
        R = hcat(R...)'
        lay = @layout [a; b{0.2h}]
        plot(heatmap(R[ind,:], clim=(0,1),c=:vik), plot(mean(R,dims=1)'),
             layout=lay)
        Plot.save((;k1,k2))
   end


   Plot.setfolder("heatmap of individual groups norm01")
   k = first(keys(keyz))
   H=[]
   for k in keys(keyz)
       key = keyz[k]
       f   = prep(matrixform(F[key]))
       shifts = vec(f.dims[2].val.data)
       inds = sortperm(f[:bestshift_bitsperspike][:,1])
       bps  = f[:bitsperspike][inds, :]
       XX = hcat([Utils.norm_extrema(b) for b in eachrow(bps)]...)'
       h=heatmap(shifts, collect(1:size(XX,1)), XX; label="",
       title="$(key.datacut)", size=(600,1200), c=:hot, clim=(0.3,1))
       vline!([0],c=:pink,linestyle=:dash, linewidth=3,label="")
       Plot.save((;k))
       push!(H,h)
   end
   plot(H...,size=(3000,3000)) 
   Plot.save("together")

   Plot.setfolder("heatmap of individual groups normP")
   k = first(keys(keyz))
   H=[]
   for k in keys(keyz)
       key = keyz[k]
       f   = prep(matrixform(F[key]))
       shifts = vec(f.dims[2].val.data)
       inds = sortperm(f[:bestshift_bitsperspike][:,1])
       bps  = f[:bitsperspike][inds, :]
       XX = hcat([Utils.norm_percent(b,0.5) for b in eachrow(bps)]...)'
       h=heatmap(shifts, collect(1:size(XX,1)), XX; label="",
        title="$(key.datacut)", size=(600,1200), c=:hot, clim=(0.3,1))
        vline!([0],c=:pink,linestyle=:dash, linewidth=3,label="")
        Plot.save((;k))
        push!(H,h)
    end
    plot(H...,size=(3000,3000)) 
    Plot.save("together")

    Plot.setfolder("heatmap of individual groups")
    k = first(keys(keyz))
    H=[]
    for k in keys(keyz)
        key = keyz[k]
        f   = prep(matrixform(F[key]))
        shifts = vec(f.dims[2].val.data)
        inds = sortperm(f[:bestshift_bitsperspike][:,1])
        bps  = f[:bitsperspike][inds, :]
        XX = hcat([b for b in eachrow(bps)]...)'
        h=heatmap(shifts, collect(1:size(XX,1)), XX; label="",
        title="$(key.datacut)", size=(600,1200), c=:hot, clim=(0.3,1))
        vline!([0],c=:pink,linestyle=:dash, linewidth=3,label="")
        Plot.save((;k))
        push!(H,h)
    end
    plot(H...,size=(3000,3000)) 
    Plot.save("together")

end # dataset


#=

statfunc = mean
nanfunc = eval(Symbol("nan"*string(Symbol(statfunc))))
Plot.setfolder("fixed_$(string(Symbol(statfunc)))maps_centered")
# time: 2022-09-08 18:05:44 EDT
# mode: julia
# Visualize the heatmap
# time: 2022-09-08 18:05:44 EDT
# mode: julia
P=[]
# time: 2022-09-08 18:05:45 EDT
# mode: julia
for (k1, k2) in zip((:nontask, :cue, :mem_error, :cue_error, :error), (:task,:memory, :mem_correct, :cue_correct, :correct))
    @info "keys" k1 k2
    key1,key2 = keyz[k1],keyz[k2]
    f1, f2 = prep(matrixform(F[key1])), prep(matrixform(F[key2]))
    u = intersect(f1.dims[1], f2.dims[1])
    f1,f2 = f1[unit=At(u)], f2[unit=At(u)]

    inds = sortperm(f1[:bestshift_bitsperspike][:,1])
    XX1 = bps  = f1[:bitsperspike][inds, :]
    XX1 = bps = XX1 .- Matrix(statfunc(XX1, dims=2))
    V1 = vec(statfunc(XX1,dims=1))
    #bs1 = bootstrap(func, collect(eachrow(XX1)), BalancedSampling(1000))
    #straps1 = hcat(straps(bs1)...)

    inds = sortperm(f2[:bestshift_bitsperspike][:,1])
    XX2 = bps  = f2[:bitsperspike][inds, :]
    XX2 = XX2 .- Matrix(statfunc(XX2, dims=2))
    V2 = vec(statfunc(XX2,dims=1))
    #bs2 = bootstrap(func, collect(eachrow(XX2)), BasicSampling(1000))
    #straps2 = hcat(straps(bs1)...)

    # I can do this because I intersected neurons in both
    #bs21 = bootstrap(func, collect(eachrow(XX2-XX1)), BasicSampling(1000))
    #ci21 = confint(bs21, BasicConfInt(0.95))
    #straps21 = [x-y for (x,y) in Iterators.product(eachrow(straps1),
                                                   #eachrow(straps2))]
    shifts = collect(f1.dims[2])
    # Manual jackknife
    M = Matrix(XX2-XX1)
    jacknives = [statfunc(M[setdiff(1:size(M,1),[i]),:],dims=1)
                 for i in 1:size(M,1)]
    plot(hcat([[y for y in x] for x in ci21]...)')
    
    p=plot(shifts, V2-V1, title="$k2 - $k1", size=0.75 .* (800,1200), fillrange=0,
          xlabel="shift", ylabel="Δ bits per spike")
    vline!([0],c=:black,linestyle=:dash, linewidth=2)
    plot!(shifts,hcat([vec(m) for m in jacknives]...),c=:black,legend=:none, linestyle=:dot, alpha=0.4,title="Jacknifed $(string(Symbol(statfunc))) difference in\n bits per spike, $k2-$k1");
    hline!([0],c=:red,linestyle=:dash)
    
    if (k1,k2) == (:cue, :memory) && statfunc == mean
        shifts_local = Utils.in_range(shifts, [-.21, .61])
        shifts_distal = (!).(shifts_local)
        get_local(val)  = nanfunc(val[shifts_local])
        get_distal(val) = nanfunc(val[shifts_distal])
        L, D = get_local((V2-V1)),
               get_distal((V2-V1))
       stde(x) = nanstd(x,dims=1)/sqrt(size(x,1))
       J = vcat(jacknives...)
       Lse, Dse = statfunc(get_local(stde(J))), statfunc(get_distal(stde(J)))
       b=bar(["local", "nonlocal"],   [L,D], yerror=[Lse, Dse], label="")
       ##scatter!(["local", "nonlocal"], [L,D], err=[Lse, Dse],
                #markerstrokecolor = :transparent,
                #linestyle=:dash, c=:black, markersize=10, label="", lw=4)
       plot!(["local", "nonlocal"], [L,D], err=[Lse, Dse],
                l=nothing, m=(4,stroke(:red, 1.5)), label="")

       
       plot(p, b, size=(1000,1000))
       Plot.save((;desc="barsum", k1, k2))
    end
    push!(P,p)
    p
    Plot.save((;xlim="full", k1, k2))
end
# time: 2022-09-08 18:05:51 EDT
# mode: julia
plot(P...)
# time: 2022-09-08 18:05:51 EDT
# mode: julia
Plot.save("together")

# Seems to uncover an effect in the mean bits per spike,
# with current being -0.2s to 0.4s


Plot.setfolder(folder,"fixed_heatmaps3")
for key in keys(F)
    f = prep(matrixform(F[key]))
    #push_metric!(f, metrics.bitsperspike)
    #push_shiftmetric!(f, best_tau!; metric=:bitsperspike)
    shifts = vec(f.dims[2].val.data)
    inds = sortperm(f[:bestshift_bitsperspike][:,1])
    bps  = f[:bitsperspike][inds, :]
    XX = hcat([Utils.norm_extrema(b) for b in eachrow(bps)]...)'
    heatmap(shifts, collect(1:size(XX,1)), XX, 
    title="$(key.datacut)", size=(600,1200), c=:hot, clim=(0.3,1))
    vline!([0],c=:black,linestyle=:dash, linewidth=2)
    Plot.save((;xlim="full", key...))
end


function field_grad(f)
    np = pyimport("numpy")
    ff = [ff.rate for ff in f]
    cell = cat(ff[2,:]...,dims=3)
    np.gradient(cell)
end

function viz_gradient(nabla, scale=1, ck=false, scalek=false)

    xx,yy,zz,ii,jj,kk = ([] for _ in 1:6)
    for (xyz, i, j, k) in zip(CartesianIndices(nabla[1]), nabla...)
        x,y,z = xyz.I
        sc = scale isa Vector ? scale[z] : scale
        push!(xx,x)
        push!(yy,y)
        push!(zz,z)
        push!(ii,i*scale)
        push!(jj,j*scale)
        push!(kk,k*scale)
    end
    plt = pyimport("matplotlib.pyplot")
    #fig,axs = plt[:subplots](Int(ceil(sqrt(maximum(zz)))), Int(ceil(sqrt(maximum(zz)))))
    
    @gif for zi in 1:160
        x,y,z,i,j, k =getindex.([xx,yy,zz,ii,jj,kk], [zz.==zi])
        #axs[zi][:quiver](x, y, i, j)
        arrw = (i, j)
        arrw = scalek ? arrw .* abs.(k) : arrw
        quiver(x,y, quiver= 0.05 .* arrw , c=:black)
    end

end
=#
