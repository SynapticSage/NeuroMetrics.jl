using ProgressMeter


# Something that should work, predicting a brain area from
# itself. It would be flipping wild if I can't do this and
# potentially indicative of a deeper problem.

XX, _ = dx_dy[(0)]
_, y  = dx_dy[(0)]
f = construct_predict_spikecount(df, cells, "CA1"; dep_area="CA1")

m = glm_mlj(f[1], XX, y, Poisson(); type=:custom)
ms = [glm_mlj(ff, XX, y, Poisson(), type=:specific) for ff in f]
mc = [glm_mlj(ff, XX, y, Poisson(), type=:custom)   for ff in f]
prog = Progress(length(f))
mm=[]
for ff in f
    push!(mm, (next!(prog); glm_matlab(ff, XX, y, Poisson())))
end
# for (ff,m) in zip(f,mm)
#     next!(prog)
#     glm_matlab(ff, XX, y, Poisson(); pre=m)
# end


#MAE
vs = filter(x->abs(x) < 10e8, getindex.(ms,["adjr2"]))
vc = filter(x->abs(x) < 10e8, getindex.(mc,["adjr2"]))
vm = filter(x->abs(x) < 10e8, getindex.(mm,["adjr2"]))
using HypothesisTests
KruskalWallisTest(vs,vc)
histogram(vs)
histogram(vc)
median(vs ) - median(vc)
median(vc)
median(vs)
median(vm)

vsy = filter(x->true, getindex.(ms,["y"]))
vcy = filter(x->true, getindex.(mc,["y"]))
vsypred = filter(x->true, getindex.(ms,["ypred"]))
vcypred = filter(x->true, getindex.(mc,["ypred"]))

P=[(plot(yy);plot!(ypred)) for (yy,ypred) in zip(vcy, vcypred)]
plot(P[1:10]...)

P=[(plot(yy);plot!(ypred)) for (yy,ypred) in zip(vsy, vsypred)]
plot(P[1:10]...)

# LIteally just one neuron predicting one neuron
xn,yn = XX[:,"8"][:,DIutils.na], XX[:,"8"]
R=  lm.PoissonRegressor()
R.fit(xn, yn)
Yn = R.predict(xn)


logval(xn) = log.(vec(replace(xn,0=>1e-16)))
# logval1(xn) = log.(vec(replace(xn,0=>1e-16)) .+ exp(0) ) # Makes it worse

using Statistics
begin
    xnz = logval(xn)
    mat"[$B,$stats] = lassoglm(double($xnz),double($yn), 'poisson', 'CV', 5)"
    idxLambdaMinDeviance = stats["IndexMinDeviance"]
    intercept = stats["Intercept"]
    B0   = intercept[Int(idxLambdaMinDeviance)];  
    c = [B0; B[ :,   Int(idxLambdaMinDeviance) ]]
    ypred = mat"glmval($c, double($(xnz)), 'log', $stats)"
    ypred2 = exp.(xnz .* B[ :,   Int(idxLambdaMinDeviance) ] .+ B0)
end
begin
    mat"pGLMwts = glmfit(double($(xn)),double($yn),'poisson', 'constant', 'on')"
    pGLMconst = pGLMwts[1];
    pGLMfilt = pGLMwts[2:end]; 
    mat"$ratepred_pGLM = exp(($pGLMconst + $xn * $pGLMfilt))"
end
p = plot(plot(ypred), plot(ypred2))
plot(yn)
MLJ.mae(vec(yn), vec(ypred))
MLJ.rsquared(vec(yn), vec(ypred))

using GLM
m =glm(xnz, yn, Poisson())
ypred = GLM.predict(m, Float64.(xnz))

plot(plot(yn), plot(ypred))

plot(plot(yn),plot(Yn),title="r2=$(adjusted_r2_score(Yn,yn,length(yn)))")


using TuringGLM

