#    _  _          _               _                     
#  _| || |_    ___| |__   ___  ___| | __ __   _____  ___ 
# |_  ..  _|  / __| '_ \ / _ \/ __| |/ / \ \ / / _ \/ __|
# |_      _| | (__| | | |  __/ (__|   <   \ V /  __/ (__ 
#   |_||_|    \___|_| |_|\___|\___|_|\_\   \_/ \___|\___|
#                                                        
# Checks theta power and velocity correlation holds
# for tetrode

if !isdefined(Main, :lfp)
    include("../load_isolated.jl")
end

beh.index = 1:size(beh,1)
DIutils.filtreg.register(beh, lfp, on="time", transfer=["index"])
lf = dropmissing(combine(groupby(lfp, :index), :amp=>mean),
                 :index)
lf.vel = abs.(beh[lf.index, :velVec])

@df @subset(lf,:vel .> 4, :vel .< 200) begin
    scatter(:vel, :amp_mean, alpha=0.02, label="", xlabel="velocity", ylabel="theta power")
end

using GLM
form =@formula amp_mean ~ 1 + vel
l = lm(form, @subset(lf, :amp_mean .< 300, :vel .> 5))
form2 =@formula amp_mean ~ 1 + vel * vel^2
l2 = lm(form2, @subset(lf, :amp_mean .< 300, :vel .> 5))

using RecipesBase
@recipe function plot(l::StatsModels.TableRegressionModel;xvar=2)
    response = l.mm.m * l.model.pp.beta0
    (l.mm.m[:,xvar], response) 
end

plot!(l, xlim=(0,80), ylim=(0,300), label=string(form))
scatter!(l2, xlim=(0,80), ylim=(0,300), label=string(form2))

