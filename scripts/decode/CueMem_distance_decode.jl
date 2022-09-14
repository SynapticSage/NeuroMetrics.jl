import Plot
using HypothesisTests
using StatsBase

doanimation =  false

ripples, cycles, beh, lfp, dat = copy.((ripples, cycles, beh, lfp, dat))
GC.gc()
@time ripples, cycles = annotate_vector_info(ripples, cycles, beh, lfp, dat, x, y, T)

Utils.filtreg.register(beh, ripples, on="time", transfer=["cuemem","period"])
Utils.filtreg.register(beh, cycles, on="time",  transfer=["cuemem","period"])

theta_length = :dec₀₁

"""
General stats
"""

# RIPPPLE LENGTHS
#violin!(fill(i, size(x)), x, c=:gray, label="")
rip = @df @subset(ripples,abs.(theta_length) .> 0) Plots.violin(:cuemem , abs.(theta_length), c=:gray)
rip = @df @subset(ripples,abs.(theta_length) .> 0) Plots.scatter!(:cuemem .+ 0.15 .* randn(size(:cuemem)), abs.(theta_length))

# THETA LENGTHS
cyc = @df @subset(cycles,abs.(theta_length) .> 0) Plots.violin(:cuemem , abs.(theta_length), c=:gray)
cyc = @df @subset(cycles,abs.(theta_length) .> 0) Plots.scatter!(:cuemem .+ 0.15 .* randn(size(:cuemem)), abs.(theta_length), alpha=0.1)
Plots.plot(rip,cyc, layout=Plots.grid(2,1), title="Ripple and theta lengths, cue vs mem")

# ECDF plot
@df ripples ecdfplot(abs.(theta_length))
@df cycles  ecdfplot!(abs.(theta_length))

var1, var2 = @subset(cycles, :cuemem .== -1), @subset(cycles, :cuemem .== 0 )
var1, var2 = abs.(var1.dec₀₁), abs.(var2.dec₀₁)
var1 = var1[var1.!=0]
var2 = var2[var2.!=0]
ecdfplot(var1)
ecdfplot!(var2)
pvalue(HypothesisTests.KruskalWallisTest(abs.(var1), abs.(var2)))

"""
CUE vs MEM
"""

var1, var2 = @subset(cycles, :cuemem .!= -1), @subset(cycles, :cuemem .== 1)
var1, var2 = abs.(var1.dec₀₁), abs.(var2.dec₀₁)
var1 = var1[var1.!=0]
var2 = var2[var2.!=0]
ecdfplot(var1)
ecdfplot!(var2)
pvalue(HypothesisTests.KruskalWallisTest(abs.(var1), abs.(var2)))

#@df @subset(cycles, :cuemem .== -1)  ecdfplot(abs.(theta_length))
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


#@df @subset(ripples, :cuemem .== -1)  ecdfplot(abs.(theta_length))
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

"""
CUE vs MEM

periodwise
"""
cycles.isofrac = Vector{Float64}(undef, size(cycles,1))
cycles.isocount = Vector{Float64}(undef, size(cycles,1))
cycles.spcount = Vector{Float64}(undef, size(cycles,1))
for cycle in eachrow(cycles)
    sp = Utils.in_range(spikes.time, [cycle.start, cycle.stop])
    cycle.isofrac, cycle.isocount = mean(spikes[sp, :isolated]), sum(spikes[sp, :isolated])
    cycle.spcount = sum(sp)
end

@df cycles Plots.scatter(:isocount, abs.(theta_length), alpha=0.2)

@df cycles Plots.histogram(:isofrac, label="")
@df cycles Plots.vline!([quantile(:isofrac, 0.99)],label="0.99")
@df cycles Plots.vline!([quantile(:isofrac, 0.95)],label="0.95")
@df cycles Plots.vline!([quantile(:isofrac, 0.90)],label="0.90")


cycles_clean = @subset(cycles, :isofrac .< 0.35)
bins = 7
cycles_clean.icbin = Utils.binning.digitize(cycles_clean.isocount, bins+1)
cycles_clean.isofrac = replace(cycles_clean.isofrac, NaN=>0)
cycles_clean.ifbin = Utils.binning.digitize(cycles_clean.isofrac, bins+1)
@subset!(cycles_clean, :ifbin .<= bins,
         abs.(theta_length) .> 0,
         abs.(theta_length) .> quantile(abs.(theta_length), 0.01),
         abs.(theta_length) .< quantile(abs.(theta_length), 0.99),
         abs.(:isofrac) .> quantile(abs.(:isofrac), 0.01),
         abs.(:isofrac) .< quantile(abs.(:isofrac), 0.99),
        )

# -------------------------------------
# BINNED PLOTS: Grouping by isofrac or isocount alone
# -------------------------------------
sumtab = combine(groupby(cycles_clean,:icbin), theta_length => (x->mean(abs.(x))) => :mean,
        theta_length => (x->median(abs.(x))) => :median
       )
@df sumtab Plots.plot(:icbin, :median)
@df sumtab Plots.plot!(:icbin, :mean)

sumtab = combine(groupby(cycles_clean,:ifbin), theta_length => (x->mean(abs.(x))) => :mean,
        theta_length => (x->median(abs.(x))) => :median
       )

@df sumtab Plots.plot(:ifbin, :median)
@df sumtab Plots.plot!(:ifbin, :mean)
@df sumtab Plots.histogram(:isofrac, title="Throw out events that are 100% iso")


# -------------------------------------
# LINEAR MODELS
# -------------------------------------
sumtab = combine(groupby(cycles_clean,[:isofrac,:spcount, :cuemem]), 
                 theta_length => (x->mean(abs.(x))) => :mean,
        theta_length => (x->median(abs.(x))) => :median
       )


"""
Full linear and quadratic model
"""
# Linear and quadratic term
using GLM
form_lin_quadratic = @formula(mean ~ isofrac * spcount * isofrac^2 * spcount)
#C = coef(l) .* [ones(size(sumtab,1)), sumtab.isofrac, sumtab.spcount, sumtab.isofrac.^2]
#l = lm(form_lin_quadratic, @subset(sumtab,:isofrac .!= maximum(sumtab.isofrac)))
l = lm(form_lin_quadratic, sumtab)
@info "r²" form_lin_quadratic adjr²(l)
Plots.scatter(modelmatrix(l)[:,2], response(l), alpha=0.1)
Plots.plot!(modelmatrix(l)[:,2], modelmatrix(l) * coef(l))

"""
R2 0.035
"""
# LINEAR - IF SPCOUNT
form_linear = @formula(mean ~ isofrac * spcount * cuemem)
l = lm(form_linear, @subset(sumtab,:isofrac .!= maximum(sumtab.isofrac)),
       contrasts=Dict(:cuemem=>DummyCoding())
      )
@info "r²" form_linear adjr²(l)

Plots.scatter(modelmatrix(l)[:,2], response(l), alpha=0.1) 
Plots.plot!(modelmatrix(l)[:,2], modelmatrix(l) * coef(l))
if doanimation
    anim = Plots.@animate for ang in 1:360
        Plots.scatter(eachcol(modelmatrix(l)[:,2:3])..., response(l), alpha=0.1, camera=(ang,15))
        Plots.plot!(eachcol(modelmatrix(l)[:,2:3])..., modelmatrix(l) * coef(l), alpha=0.2)
    end
    Plots.gif(anim)
end


"""
R2 0.005
"""
# LINEAR - IF
form_isofrac = @formula(mean ~ isofrac)
l = lm(form_isofrac, @subset(sumtab,:isofrac .!= maximum(sumtab.isofrac)))
apply_schema(form_isofrac, schema(sumtab))
Plots.scatter(modelmatrix(l)[:,2], response(l), alpha=0.1)
Plots.plot!(modelmatrix(l)[:,2], modelmatrix(l) * coef(l), xlabel="fraction isolated", ylabel="θ distance")
@info "r²" form_isofrac adjr²(l)



"""
R2 0.0153
"""
# LINEAR - IF
replace(sumtab.cuemem, missing=>-2)
form_isofrac = @formula(mean ~ isofrac + cuemem)
l = lm(form_isofrac, @subset(sumtab,:isofrac .!= maximum(sumtab.isofrac)), contrasts=Dict(:cuemem=>DummyCoding()))
apply_schema(form_isofrac, schema(sumtab))
c = get(ColorSchemes.vik, Utils.norm_extrema(disallowmissing(sumtab.cuemem)))
Plots.scatter(modelmatrix(l)[:,2], response(l), alpha=0.1, c=c)
Plots.plot!(modelmatrix(l)[:,2], modelmatrix(l) * coef(l), xlabel="fraction isolated", ylabel="θ distance")
@info "r²" form_isofrac adjr²(l)
