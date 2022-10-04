
using Distributed
using ProgressMeter
addprocs(2)
#addprocs([("mingxin",2),])
@everywhere using DrWatson
@everywhere quickactivate(expanduser("~/Projects/goal-code"))
@everywhere using PyCall
using ThreadSafeDicts
using DataFramesMeta
pids = workers()

# Load data
# ----------------
@time spikes, beh, ripples, cells = Load.load("RY16", 36);

@everywhere areas=(:ca1,:pfc)
R  = Dict()
R[:ca1], R[:pfc] = (Munge.spiking.torate(@subset(spikes,:area .== ar), beh)
                    for ar in ("CA1","PFC"))


method = :FisherS
method = :TwoNN


# GLOBAL intrinsice
#prog = Progress(2;desc="global fit")
G =  []
for area in areas
    g= begin
        # Intrinsic dimension
        skd = pyimport("skdim")
        dimMeasure = getproperty(skd.id, method)()
        dimMeasure.fit(Matrix(R[area] .+ eta)).dimension_
    end
    push!(G,g)
end
globaldim = Dict(a=>g for (a,g) in zip(areas, G))

# POINTWISE local intrinsic
L = []
#prog = Progress(2;desc="local fit")
for area in areas
    l = begin
        # Intrinsic dimension
        skd = pyimport("skdim")
        #dimMeasure = skd.id.lPCA()
        dimMeasure = skd.id.FisherS()
        dimMeasure.fit_pw(R[area], n_neighbors=50, n_jobs=1).dimension_pw_;
    end
    push!(L,l)
end
localdim = Dict(a=>l for (a,l) in zip(areas, L))

using Serialization
serialize(datadir("manifold", "intrinsic_method=$method.serial"),
          (;localdim,globaldim))

D = DataFrame([localdim[:pfc] localdim[:ca1] beh.cuemem],
              [:ld_pfc, :ld_ca1, :cuemem])


globaldim, localdim = Dict(k=>(try;fetch(v);catch;v;end) for (k,v) in globaldim),
                      Dict(k=>(try;fetch(v);catch;v;end) for (k,v) in localdim)
                      
beh.dimca1, beh.dimpfc = localdim[:ca1], localdim[:pfc]
numcellsca1,numcellspfc = sum(R[:ca1] .> 0, dims=2),
                          sum(R[:pfc] .> 0, dims=2)

p1 = plot(numcellspfc, beh.dimpfc, alpha=0.01, xlabel="neurons", ylabel="dims")
plot!(1:15,1:15,lw=2,c=:black,linestyle=:dot, title="pfc\n globaldim=$(round(globaldim[:pfc],sigdigits=2))")
p2 = plot(numcellsca1, beh.dimca1, alpha=0.01, xlabel="neurons")
plot!(1:40,1:40,lw=2,c=:black,linestyle=:dot, title="ca1\n globaldim=$(round(globaldim[:ca1],sigdigits=2))")
plot(p1,p2,aspectratio=1)

savefig(plotsdir("manifold", "intrinsic_dimension_area_local_method=$method.png"))
savefig(plotsdir("manifold", "intrinsic_dimension_area_local_method=$method.pdf"))


D = DataFrame([localdim[:pfc] localdim[:ca1] numcellspfc numcellsca1 beh.cuemem],
              [:ld_pfc, :ld_ca1, :nc_ca1, :nc_pfc, :cuemem])

combine(groupby(D, :cuemem), :ld_ca1=>x->nanmean(collect(Utils.skipnan(x))))
combine(groupby(D, :cuemem), :ld_pfc=>nanmean)
combine(groupby(D, :cuemem), :ld_pfc=>x->nanstd(x)/sqrt(length(x)))
combine(groupby(D, :cuemem), :ld_ca1=>x->nanstd(collect(Utils.skipnan(x)))/sqrt(length(x)))

combine(groupby(D, :cuemem), :nc_ca1=>x->nanmean(collect(Utils.skipnan(x))))
combine(groupby(D, :cuemem), :nc_pfc=>nanmean)
combine(groupby(D, :cuemem), :nc_pfc=>x->nanstd(x)/sqrt(length(x)))
combine(groupby(D, :cuemem), :nc_ca1=>x->nanstd(collect(Utils.skipnan(x)))/sqrt(length(x)))
