

# Intrinsic dimension
skd = pyimport("skdim")
#dimMeasure = skd.id.lPCA()
dimMeasure = skd.id.FisherS()
# GLOBAL intrinsice
globaldim = Dict()
prog = Progress(2;desc="global fit")
for area in areas
    globaldim[area] = dimMeasure.fit(R[area]).dimension_
    next!(prog)
end
# POINTWISE local intrinsic
localdim = Dict()
prog = Progress(2;desc="local fit")
for area in areas
    localdim[area] = dimMeasure.fit_pw(R[area], n_neighbors=50,
                                       n_jobs=8).dimension_pw_
    next!(prog)
end
globaldim, localdim = Dict(k=>fetch(v) for (k,v) in globaldim),
                      Dict(k=>fetch(v) for (k,v) in localdim)
beh.dimca1, beh.dimpfc = localdim[:ca1], localdim[:pfc]
numcellspfc,numcellsca1 = sum(R[:ca1] .> 0, dims=2),
                          sum(R[:pfc] .> 0, dims=2)

p1 = plot(numcellspfc, beh.dimpfc, alpha=0.01, xlabel="neurons", ylabel="dims")
plot!(1:15,1:15,lw=2,c=:black,linestyle=:dot)
p2 = plot(numcellsca1, beh.dimca1, alpha=0.01, xlabel="neurons")
plot!(1:40,1:40,lw=2,c=:black,linestyle=:dot)
plot(p1,p2,aspectratio=1)
savefig(plotsdir("manifold", "intrinsic_dimension_area_local.png"))
savefig(plotsdir("manifold", "intrinsic_dimension_area_local.pdf"))

