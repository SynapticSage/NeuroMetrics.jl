using Interact, Blink, Mux, ProgressMeter, Statistics, DataStructures
using DataFrames
animal, day, abbreviated = "RY16", 36, true
if abbreviated
    data_source=[ "spikes","behavior"]
    spikes, beh  = raw.load(animal, day; data_source=data_source)
else
    data_source=[ "spikes","behavior","lfp","cells","tetrode","ripples"]
    spikes, beh, lfp, cells, tetrode, ripples = raw.load(animal, day; 
                                                        data_source=data_source)
end

# Filter
bprops = ["poke", "correct", "velVec"];
filters = Dict("velVec"=>x->abs.(x) .≤ 0.5)
lookupcols = Dict((source=1,target=2) =>
                  table.expand_colnames(beh,bprops))
beh, spikes = raw.filterTables(beh, spikes; filters=filters,
                                            lookupcols=lookupcols)
spikes[!,"poke"] = spikes[!,"poke_1"] .+ spikes[!,"poke_4"] .+ 
                   spikes[!,"poke_2"] .+ spikes[!,"poke_5"] .+
                   spikes[!,"poke_3"]
r(x) = replace(x,NaN=>0)
beh[!,"poke"] = r(beh[!,"poke_1"]) .+ r(beh[!,"poke_4"]) .+ 
        r(beh[!,"poke_2"]) .+ r(beh[!,"poke_5"]) .+
        r(beh[!,"poke_3"])
S = spikes[spikes[!,"poke"].>=1,:]
B = beh[beh[!,"poke"].>=1,:]
B = combine(groupby(B,"correct"), nrow=>:count)
B_ratio = B[B.correct.==1,:].count ./ B[B.correct.==0,:].count

# Find best reward
units = groupby(S, ["unit"])
D = DataFrame()
for unit in units
    G = combine(groupby(unit, ["unit","correct"]), nrow=>:count)
    G = G[G.correct.>=0,:]
    if sum(G.correct .== 0) == 0
        append!(G, DataFrame(G[end,:]))
        G[end,"correct"] = 0
        G[end,"count"] = 0
    end
    if sum(G.correct .== 1) == 0
        append!(G, DataFrame(G[end,:]))
        G[end,"correct"] = 1
        G[end,"count"] = 0
    end
    G[:,"ratio"] .= G.count[G.correct.==1,:]./G.count[G.correct.==0,:]
    G[:,"ratio_norm"] = G[:,"ratio"] ./ B_ratio
    append!(D,G)
end

D = D[(D.correct .== 1) .& (D.count .> 100), :]
D = sort(D, "ratio_norm", rev=true)

# PSTHs
pokeTimes = beh[findall(diff(beh.poke).==1),:]
pokeTimes = pokeTimes[pokeTimes.correct.==1,:]
pokeTime = pokeTimes.time[1]
units = D[1:30,:unit]
n_boot = 1000

function plot_cell(unit, spikes, pokeTimes)

    start = (pokeTimes.time.-0.5)
    stop  = (pokeTimes.time.+1.5)
    timegrid = [range(start[i], stop[i], length=40) for i in 1:length(start)]
    mtimegrid = mean(hcat(collect.(timegrid)...)' .- pokeTimes.time ,dims=1)
    mtimegrid = mean([mtimegrid[1:end-1];; mtimegrid[2:end]],dims=2)
    Δ = mtimegrid[2]-mtimegrid[1]
    hists = fit.(Histogram, [spikes[spikes.unit.==unit,:].time], collect.(timegrid))
    hists = hcat([hist.weights for hist in hists]...)
    filt = (!).(all(hists.==0, dims=1))
    hists_filt = hists[:,vec(filt)]
    
    tot(x) = vec(mean(x,dims=1))
    bs1 = Bootstrap.bootstrap(tot, collect(hists'), BasicSampling(n_boot))
    cil = 0.95;
    bci1 = confint(bs1, BasicConfInt(cil));
    sums = vcat([hcat(x...) for x in bci1]...)
    ribbon = (abs.(sums[:,2]-sums[:,1]), abs.(sums[:,3]-sums[:,1]))

    return Plots.plot(
               Plots.heatmap(vec(mtimegrid), vec(1:size(hists_filt,2)),
                             hists_filt', xlabel="time", ylabel="trial"),
               Plots.plot(vec(mtimegrid), sums[:,1], label="Total spikes",
                          xlabel="time", ribbon=ribbon), layout=grid(2,1)

              )
end


w = Window();
ui = @manipulate for unit in units
    plot_cell(unit, spikes, pokeTimes)
end
body!(w,ui)

for unit in [18, 45, 46, 48, 50, 54]
    p = plot_cell(unit, spikes, pokeTimes)
    savefig(p,plotsdir("reward_responsive", "cell_$unit.svg"))
    savefig(p,plotsdir("reward_responsive", "cell_$unit.png"))
end

maxFr, meanFr, medianFr = [],[],[]
for cell in groupby(spikes, "unit")
    mintime = minimum(cell.time)
    maxtime = maximum(cell.time)
    G = range(mintime, maxtime, length=1000)
    Δ = G[2]-G[1]
    X = fit(Histogram, cell.time, collect(G)).weights
    x = maximum(X)./Δ
    y = mean(X)./Δ
    z = median(X)./Δ
    push!(maxFr, x)
    push!(meanFr, y)
    push!(medianFr, z)
end
Plots.histogram(meanFr, label="Mean Firing Rates", xscale=:identity, xticks=xlims())
Plots.histogram(meanFr, label="Mean Firing Rates", xscale=:identity, xticks=xlims())
