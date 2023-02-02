@time spikes, beh, ripples, cells = Load.load("RY16", 36);
Load.filtreg.register(beh,spikes,on="time", transfer=["traj"])
beh = Load.filtreg.filter(beh, filters=filt[:all], filter_skipmissingcols=true)
spikes = Load.filtreg.filter(spikes, filters=filt[:all], filter_skipmissingcols=true)

B = Dict(g.traj[1] =>
     (p=plot(g.time, g.x); plot!(g.time,g.y); p)
     for g in groupby(beh, :traj))
B[207]

#using Munge.SpikeTrains

S = Dict(s.traj[1] => (p=scatter(s.time,s.unit,m=:vline))
         for  s in groupby(spikes,:traj))
S[207]
lay = @layout [a; b{0.2h}]

for ind in keys(S)
    if  ismissing(ind)
        continue
    end
    plot(S[ind], B[ind], layout=lay)
    savefig("/home/ryoung/Projects/goal-code/plots/behavior/trajectBeh-with-spiking/xy-$ind.pdf")
end

# fake behavior
using Distributions
G=  Distributions.Gaussian(0, 1)
    g = 2 .* pdf.([G], collect(-3:0.1:3))

for i in 100:100:(1000-length(g))
    Bf = rand(10, 1000)
    start, final = i, min((i + length(g)-1), size(Bf,2))
    vBf =  view(Bf, :, start:final) 
    vBf .= vBf .+ g[Utils.na, :]
    plot(
         heatmap(Bf, colorbar=false),
        plot(Bf', xlabel="time", ylabel="trial", alpha=0.6, legend=:none),
        layout=grid(2,1)
       )
    savefig(plotsdir("behavior","fake-behavior","$(i/1000) perecent.pdf"))
end
