using Serialization
using ProgressMeter
import Plot
using DataStructures: OrderedDict
using Plots
using HypothesisTests
using SoftGlobalScope
using DataFrames
using Timeshift
using Timeshift.shiftmetrics
using Field.metrics
using Statistics, NaNStatistics, HypothesisTests
using DimensionalData
using JLD2

se(x)  = std(x)/sqrt(length(x))
datacut = :all
w = 8
animal,day="RY16",36
#f,fwg = deserialize(datadir("exp_pro", "xyG-$datacut-$w-fmat"));
f, fwg, animal, day, shifts, thresh, w = 
    deserialize(datadir("exp_pro", "$animal-$day-xyG-$datacut-$w-fmat"), )


# Add metrics
push_shiftmetric!(fwg, best_tau!; metric=:bitsperspike)
push_dims!(fwg)
push_dims!(f)

inds = vec(any(fwg[:meanrate] .> 0.001 .&& 
               fwg[:area] .== "CA1",dims=2))

f, fwg = f[inds, :], fwg[inds,:]


# Get shifts
sh = collect(shifts)
# -------------
# GOAL INDEXING
# -------------
# Grab the goalindex as the mean response to goal-pursuit across the
# different XY positions 
gd = [Utils.squeeze(nanmean(nanmean(f.rate; dims=1), dims=2))
      for f in fwg]
goal_index_unocc = [nanmaximum(g)-nanminimum(g) for g in gd]
g_count = [Utils.squeeze(nansum(nansum(f.count; dims=1), dims=2))
            for f in fwg]
g_occ_count = [Utils.squeeze(nansum(nansum(f.occ.count; dims=1), dims=2))
            for f in fwg]
#goal_index = [nanmaximum(c./o)-nanminimum(c./o)
#      for (c,o) in zip(g_count, g_occ_count)]
goal_index = [
              (x=nanmaximum(c./o) - nanminimum(c./o);
               if sum(o).>=2 && !isinf(abs(x))
                   x
               else
                   NaN
               end
              )
              for (c,o) in zip(g_count, g_occ_count)
             ]
goals = [
              c./o 
              for (c,o) in zip(g_count, g_occ_count)
             ]

# CLEAN :  Throw away responses with best BPS less than 0.5 (papers throw out these neurons)
fwg[:bestshift_bitsperspike][fwg[:bestshift_bitsperspike] .< 0.5] .= NaN


# Obtain a dimarray representation
goal_index_unocc = DimArray(goal_index, fwg.dims)
goal_index= DimArray(goal_index, fwg.dims)

# Create a table mapping shifts to gi
shift_vs_gi = []
for (r_fwg, r_gm) in zip(eachrow(fwg),eachrow(goal_index))
    loc = r_fwg[:bestshift_bitsperspike][1]
    samp = [r_fwg[shift=At(loc)][:bestshift_bitsperspike],
             r_gm[shift=At(loc)]]
    push!(shift_vs_gi, samp)
end
shift_vs_gi = hcat(shift_vs_gi...)'
# Throw out outliers
outlier = quantile(shift_vs_gi[:,2], 0.98)
nonoutlier = shift_vs_gi[:,2] .< outlier .&& (!).(isinf.(shift_vs_gi[:,1]))
shift_vs_gi = shift_vs_gi[nonoutlier, :]
# Create sym and asym versions
shift_vs_gi_asym = copy(shift_vs_gi)
shift_vs_gi[:,1] = abs.(shift_vs_gi[:,1]) # absolute distance of τ
shift_vs_gi_asym = shift_vs_gi_asym[abs.(shift_vs_gi_asym[:,1]) .> 0,:]

sgb = []
for (i,(r_fwg, r_gm)) in enumerate(zip(eachrow(fwg),eachrow(goal_index)))
    samp = hcat(vec.([i*ones(size(r_gm)), r_fwg[:shift], r_gm, r_fwg[:bitsperspike]])...)
    push!(sgb, samp)
end
sgb = vcat(sgb...)

dogif = false
if dogif
c = Utils.plotutils.getplotcolor(sgb[:,end-1],:vik)
anim = @animate for theta in 1:2:360
    p = plot([],[],[], xlabel="τ shift\n(s)", ylabel="goal\nindex", zlabel="bits/spike", camera=(theta,30))
    for i in unique(sgb[:,1])
        inds = sgb[:,1] .== i
        s = sgb[inds,2:end]
        ci = c[inds]
        plot!(p,eachcol(s)...; c=ci, markersize=3,legend=nothing)
    end
    p
end
gif(anim, datadir(Plot.folder_args..., "shift-goalind-bps.gif"))
#Plot.setfolder("")
end


#using PyCall
#plt = pyimport("matplotlib.pyplot")
#@pyimport numpy
#np = numpy
#plt.ion()
#fig, ax = plt.subplots(subplot_kw=Dict(:projection=>"3d"))
#x,y,z = np.atleast_2d.(np.array.(vec.(collect(eachcol(sgb)))))
#ax.plot_surface(x,y,z)


partition = shift_vs_gi[:,1] .> 0.5
chunks = [(shift_vs_gi[(!).(partition),2]), (shift_vs_gi[partition,2])]
t = OneWayANOVATest(chunks...)
bar(1:2, mean.(chunks), yerror=se.(chunks), title="pval=$(pvalue(t))")

#Plot.save("summary, datacut=$datacut, w=$w, pval=$(pvalue(t))")

partition = shift_vs_gi[:,1] .> 1
chunks = [(shift_vs_gi[(!).(partition),2]), (shift_vs_gi[partition,2])]
t = OneWayANOVATest(chunks...)
bar(1:2, mean.(chunks), yerror=se.(chunks), title="pval=$(pvalue(t))")

#Plot.save("summary, <1, datacut=$datacut, w=$w, pval=$(pvalue(t))")
##scatter(vec(fwg[:bestshift_bitsperspike]), vec(gm), title="$datacut")

#histogram2d(vec(fwg[:bestshift_bitsperspike]), vec(goal_index_unocc), title="$datacut", xlabel="time", ylabel="goal-index\n\nmax FR\$_g\$-min FR\$_g\$", colorbar=nothing)
using GLM

# SCATTER NON-ZERO TIMES
type = :sym
if type == :asym
    xval = vec(shift_vs_gi_asym[:,1])
    yval = vec(shift_vs_gi_asym[:,2])
else
    xval = vec(shift_vs_gi[:,1])
    yval = vec(shift_vs_gi[:,2])
end
m = minimum(abs.(xval)) 
xval, yval = xval[abs.(xval) .> m .&& yval .< 0.5], yval[abs.(xval) .> m .&& yval .< 0.5]
X = DataFrame([xval yval], [:shift, :gi])
fitobj = fit(LinearModel, @formula(gi ~ shift), X)
if type==:asym
    xval_plus_gaus = clamp.(xval .+ 0.02 .* randn(length(xval)), -100, 10000.0)
else
    xval_plus_gaus = clamp.(xval .+ 0.02 .* randn(length(xval)), 0.002, 10000.0)
end
scatobj = scatter(xval_plus_gaus, yval, title="$datacut", markersize=4, c=:skyblue,label="")
C=coeftable(fitobj)
#plot!(xval_plus_gaus, C.cols[1][1] .+ xval_plus_gaus*C.cols[1][2]; c=:black, linestyle=:dot, alpha=0.4)
vline!([0.25],c=:black, linestyle=:dot, linewidth=2, label="")
scatobjlog = plot(scatobj, xscale=:log)

partition = xval .> 0.25 # where is our boundary
chunks = [(yval[(!).(partition)]), (yval[partition])]
t = OneWayANOVATest(chunks...)
b25 = bar(1:2, mean.(chunks), yerror=se.(chunks), title="pval=$(pvalue(t))")

Plot.save("summary <0,25, datacut=$datacut, w=$w, pval=$(pvalue(t))")

#Plot.save((;desc="scatter goal_index versus bestshift_bitsperspike",
           #datacut))

histogram(vec(fwg[:bestshift_bitsperspike]), title="$datacut", normalize=:probability, bar_width=0.07, linecolor=nothing)
vline!([0],c=:black)
vline!([-0.25, 0.25],c=:gray)
Plot.save((;desc="histogram bestshift_bitsperspike", datacut))

inds = sortperm(fwg[:bestshift_bitsperspike][:,1])
bps  = fwg[:bitsperspike][inds, :]

partition = shift_vs_gi_asym[:,1] .> 0 # where is our boundary
chunks = [(shift_vs_gi_asym[(!).(partition),2]), (shift_vs_gi_asym[partition,2])]
t = OneWayANOVATest(chunks...)
bar(["future","past"], mean.(chunks), yerror=se.(chunks), title="pval=$(pvalue(t))")
#Plot.save("future vs past, datacut=$datacut, w=$w, pval=$(pvalue(t))")

plot(scatobj, scatobjlog, b25)

#X = hcat([Utils.norm_extrema(b) for b in eachrow(bps)]...)'
#heatmap(sh, collect(1:size(X,1)), X, title="$datacut", size=(400,1000))
#vline!([0],c=:black,linestyle=:dash, linewidth=2)
#Plot.save((;desc="snake plot, norm 01", datacut))

#X = hcat([Utils.norm_percent(b,0.5) for b in eachrow(bps)]...)'
#heatmap(sh, collect(1:size(X,1)), X, title="$datacut", size=(400,1000))
#vline!([0],c=:black,linestyle=:dash, linewidth=2)
#Plot.save((;desc="snake plot, norm percent", datacut))
