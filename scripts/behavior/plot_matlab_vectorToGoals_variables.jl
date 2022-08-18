using DrWatson
quickactivate(expanduser("~/Projects/goal-code/"))
using LazyGrids, MAT
using Plots

cd(datadir("exp_pro", "behavior", "RY16-36-02_behaviorVars.mat"));
q=MAT.matopen("test.mat"); 
gV = read(q, "goalVec");
(X,Y) = ndgrid(x,y)
x=1:size(gV,1); y=1:size(gV,2)
c = maximum(abs.(gV[:,1:5]))
plot(); 
quiver(c*vec(X[:,1:5]),c*vec(Y[:,1:5]), 
       quiver=(vec(real(gV[:,1:5])), vec(imag(gV[:,1:5]))), alpha=0.15, 
       arrow=arrow(:closed,:head,0.00001,0.00001));xlims!(0*c,2000*c)
xlim(c*0, c*2000);
savefig(plotsdir("behavior", "goal-vector", "ry16_36_2_vector_sample.png"))
