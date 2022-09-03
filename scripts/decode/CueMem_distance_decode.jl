using HypothesisTests

@time ripples, cycles = annotate_vector_info(ripples, cycles, beh, lfp, dat, x, y, T)

Utils.filtreg.register(beh, ripples, on="time", transfer=["cuemem"])
Utils.filtreg.register(beh, cycles, on="time", transfer=["cuemem"])

#violin!(fill(i, size(x)), x, c=:gray, label="")
rip = @df @subset(ripples,abs.(:dec₀₁) .> 0) Plots.violin(:cuemem , abs.(:dec₀₁), c=:gray)
rip = @df @subset(ripples,abs.(:dec₀₁) .> 0) Plots.scatter!(:cuemem .+ 0.15 .* randn(size(:cuemem)), abs.(:dec₀₁))
cyc = @df @subset(cycles,abs.(:dec₀₁) .> 0) Plots.violin(:cuemem , abs.(:dec₀₁), c=:gray)
cyc = @df @subset(cycles,abs.(:dec₀₁) .> 0) Plots.scatter!(:cuemem .+ 0.15 .* randn(size(:cuemem)), abs.(:dec₀₁), alpha=0.1)
Plots.plot(rip,cyc)

@df ripples ecdfplot(abs.(:dec₀₁))
@df cycles  ecdfplot!(abs.(:dec₀₁))

var1, var2 = @subset(cycles, :cuemem .== -1), @subset(cycles, :cuemem .== 0 )
var1, var2 = abs.(var1.dec₀₁), abs.(var2.dec₀₁)
var1 = var1[var1.!=0]
var2 = var2[var2.!=0]
ecdfplot(var1)
ecdfplot!(var2)
pvalue(HypothesisTests.KruskalWallisTest(abs.(var1), abs.(var2)))

var1, var2 = @subset(cycles, :cuemem .!= -1), @subset(cycles, :cuemem .== 1)
var1, var2 = abs.(var1.dec₀₁), abs.(var2.dec₀₁)
var1 = var1[var1.!=0]
var2 = var2[var2.!=0]
ecdfplot(var1)
ecdfplot!(var2)
pvalue(HypothesisTests.KruskalWallisTest(abs.(var1), abs.(var2)))

#@df @subset(cycles, :cuemem .== -1)  ecdfplot(abs.(:dec₀₁))
var1, var2 = @subset(cycles, :cuemem .== -1), @subset(cycles, :cuemem .!= -1)
var1, var2 = abs.(var1.dec₀₁), abs.(var2.dec₀₁)
var1 = var1[var1.!=0]
var2 = var2[var2.!=0]
ecdfplot(var1)
ecdfplot!(var2)
pvalue(HypothesisTests.KruskalWallisTest(abs.(var1), abs.(var2)))

var1, var2 = @subset(cycles, :cuemem .== 0), @subset(cycles, :cuemem .== 1)
var1, var2 = abs.(var1.dec₀₁), abs.(var2.dec₀₁)
var1 = var1[var1.!=0]
var2 = var2[var2.!=0]
ecdfplot(var1)
ecdfplot!(var2)
pvalue(HypothesisTests.KruskalWallisTest(abs.(var1), abs.(var2)))

var1, var2, var3 = @subset(cycles, :cuemem .!= -1), @subset(cycles, :cuemem .== 0), @subset(cycles, :cuemem .== 1)
var1, var2, var3 = abs.(var1.dec₀₁), abs.(var2.dec₀₁), abs.(var3.dec₀₁)
var1 = var1[var1 .!= 0]
var2 = var2[var2 .!= 0]
var3 = var3[var3 .!= 0]
Makie.Figure(title="theta sequence length")
Plots.plot(ecdf(var1), label="nontask", title="theta sequence length")
Plots.plot!(ecdf(var2), label="cue")
Plots.plot!(ecdf(var3), label="memory")
Plot.setfolder("nonlocality", "thetaseq-length")
Plot.save("differential theta sequence length")


#@df @subset(ripples, :cuemem .== -1)  ecdfplot(abs.(:dec₀₁))
var1, var2 = @subset(ripples, :cuemem .== -1), @subset(ripples, :cuemem .!= -1)
var1, var2 = abs.(var1.dec₀₁), abs.(var2.dec₀₁)
var1 = var1[var1.!=0]
var2 = var2[var2.!=0]
ecdfplot(var1)
ecdfplot!(var2)
pvalue(HypothesisTests.KruskalWallisTest(abs.(var1), abs.(var2)))

var1, var2 = @subset(ripples, :cuemem .== 0), @subset(ripples, :cuemem .== 1)
var1, var2 = abs.(var1.dec₀₁), abs.(var2.dec₀₁)
var1 = var1[var1.!=0]
var2 = var2[var2.!=0]
ecdfplot(var1)
ecdfplot!(var2)
pvalue(HypothesisTests.KruskalWallisTest(abs.(var1), abs.(var2)))
