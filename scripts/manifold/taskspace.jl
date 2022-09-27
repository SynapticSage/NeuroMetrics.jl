
#  ================
# Task space projoections
# (not so great -- try mante and susillo?)
#  ================

embeddingH = Dict()
areas = (:ca1,:pfc)
embeddingH[:ca1] = umap((Rca1')[:,1:nsamp], 8)
embeddingH[:pfc] = umap((Rpfc')[:,1:nsamp], 8)

proj(x,y) = x * y' * (y * y') * y

using OneHotArrays

task_vars = replace([beh.cuemem beh.correct beh.stopWell],NaN=>-1)[1:nsamp,:]
ca1_proj = []
for i in 1:1000:(size(embeddingH[:ca1],2)-1000)
    push!(ca1_proj, 
          proj(embeddingH[:ca1][:,i:(i+1000-1)], task_vars[i:(i+1000-1),:]')
         )
end
ca1_proj = hcat(ca1_proj...)
ca1_proj = umap(ca1_proj, 3)
scatter(eachrow(ca1_proj)...)

# ---------------------------
using GLM
lm(@formula())
dfca1 = DataFrame(Rca1)
dfca1 = unstack(dfca1, :time, :unit, "", allowduplicates=true)
Utils.filtreg.register(beh, dfca1, on="time",transfer=["cuemem", "correct","stopWell"])

# NEED  TO FILL IN THE BETAS ##


##  TEMPORARY ##
form = Term(:cuemem) ~ sum(term.([x for x in propertynames(dfca1) if x ∉ [:time,:cuemem,:correct]]))
linear_cuemem = lm(form, dfca1)
adjr2(linear_cuemem)
lmcm = DataFrame(coeftable(linear_cuemem))
sum(lmcm[:,"Pr(>|t|)"] .< (0.05/size(lmcm,1)))


form = Term(:correct) ~ sum(term.([x for x in propertynames(dfca1) if x ∉ [:time,:cuemem,:correct]]))
linear_correct = lm(form, dfca1)
adjr2(linear_correct)
lmcm = DataFrame(coeftable(linear_correct))
sum(lmcm[:,"Pr(>|t|)"] .< (0.05/size(lmcm,1)))
##  TEMPORARY ##

