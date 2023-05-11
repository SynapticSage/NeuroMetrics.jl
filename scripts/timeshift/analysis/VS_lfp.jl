#=
===================================================
Script Name     : VS_lfp.jl
Author          : Ryan Young
Created         : 05-03-2023
Last Modified   : 05-03-2023
===================================================

# Description
This script plots local field potential properties against timeshifted
optimal 𝛕 of each cell.
=#

using GoalFetchAnalysis
using DataFrames, DataFramesMeta, ProgressMeter, StatsBase, Statistics,
StatsPlots, Plots, LaTeXStrings
using LaTeXStrings
import Random

#    _  _     _                    _       _       _        
#  _| || |_  | |    ___   __ _  __| |   __| | __ _| |_ __ _ 
# |_  ..  _| | |   / _ \ / _` |/ _` |  / _` |/ _` | __/ _` |
# |_      _| | |__| (_) | (_| | (_| | | (_| | (_| | || (_| |
#   |_||_|   |_____\___/ \__,_|\__,_|  \__,_|\__,_|\__\__,_|
#                                                           
animal, day = "super_clean", 0
CELLS  = DI.load_cells(animal, day, "*")
time_factors = DI.behavior_0time_before_process("super")
taus = names(CELLS)[occursin.("1of20", names(CELLS))]

# -------------------------------
# Copy pasta from setup_pyr_iso.jl
# -------------------------------
function get_data(animal, day)
    # l_pyr = DI.load_lfp(animal, day; append="pyr")
    # l_pyr = transform(l_pyr, :time=> t-> t .- time_factors[animal], renamecols=false)
    # transform!(l_pyr, :time=> t-> t .+ time_factors[animal], renamecols=false)
    cycles = DI.load_cycles(animal, day, "pyr")
    cycles.start = cycles.start_function
    cycles.stop  = cycles.stop_function
    cycles = transform(cycles,
        :start=> t-> t .- time_factors[animal],
        :stop=> t-> t .-  time_factors[animal],
        renamecols=false
    )
    cycles = cycles[!,Not([:start_function,:stop_function])]
    spikes = transform(DI.load_spikes(animal, day, "pyr_cycles_isolated"),
        :time=> t-> t .- time_factors[animal],
        renamecols=false
    )
    # spikes = subset(SPIKES, :animal=>a->a.==animal, :day=>d->d.==day)
    # Print out extrema(time) for each of these just to be sure
    begin
        println("spikes.time: ", extrema(spikes.time))
        # println("lfp.time: ", extrema(l_pyr.time))
        println("cycles.time: ", extrema(cycles.start))
    end
    DIutils.pushover("Loaded checkpoint for $animal")
    (;cycles, spikes)
end

CYCLES, SPIKES = [], [], []
datasets = [("RY16",36), ("RY22",21)]
(animal, day) = first(datasets)
for (animal, day) in datasets
    if !isfile(DI.cyclepath(animal, day, "pyr"))
        continue
    end
    println("Loading $animal")
    cycles, spikes = get_data(animal, day)
    println("Adding $taus to $animal")
    spikes = sort!(spikes, [:unit, :time])
    DIutils.filtreg.register(CELLS, spikes, on="unit", transfer=taus)
    println("Done with $animal")
    push!(CYCLES, cycles)
    push!(SPIKES, spikes)
    # Transfer the tau values to the spikes
end
CYCLES = vcat(CYCLES...)
SPIKES = vcat(SPIKES...)

#    _  _        _                _               
#  _| || |_     / \   _ __   __ _| |_   _ _______ 
# |_  ..  _|   / _ \ | '_ \ / _` | | | | |_  / _ \
# |_      _|  / ___ \| | | | (_| | | |_| |/ /  __/
#   |_||_|   /_/   \_\_| |_|\__,_|_|\__, /___\___|
#                                   |___/         
#
# INFO: Some ideas to try
# - rank instead of actual
# - individual spikes instead of mean of their phases
# - allowing different types of tau to co-occur
# - 

# First, let's actually create views into the spikes that actually occur
# in ripple and valid theta cycle
# ISSUE: better filter theta cycles before this step?
# -e.g. why theta in ripples -- violates assumption 
# - restrict to cells with significant coupling to the theta rhythm
# - why some theta phases WAY over 1 or way less than 0

SPIKES[replace(SPIKES.theta .== 0, missing=>false), :theta_phase] .= missing
SPIKES[!,:ripple] = convert(Vector{Union{Missing,Int32}}, SPIKES[!,:ripple])
SPIKES[!,:ripple_phase] = convert(Vector{Union{Missing,Float32}}, SPIKES[!,:ripple_phase])
SPIKES[replace(SPIKES.ripple .== 0,missing=>false), :ripple_phase] .= missing

rmis(x) = replace(x, missing=>false)
badinds = rmis(SPIKES.theta_phase .> 1) .|| rmis(SPIKES.theta_phase .< 0) .||
rmis(SPIKES.ripple_phase .> 1) .|| rmis(SPIKES.ripple_phase .< 0)
SPIKES = SPIKES[.!badinds,:]


SPIKES.tp = 2π*SPIKES.theta_phase
SPIKES.rp = 2π*SPIKES.ripple_phase
plot(histogram(SPIKES.tp|>skipmissing|>collect, bins=100, normed=true, label="Theta phase"),
     histogram(SPIKES.rp|>skipmissing|>collect, bins=100, normed=true, label="Ripple phase"))
SPIKES.tp = exp.(im*SPIKES.tp)
SPIKES.rp = exp.(im*SPIKES.rp)

function plot_spike_phases(sp::DataFrame)
    r = Random.randperm(nrow(sp)) 
    n=min(1000, nrow(sp))
    samples = sp[r[1:n],:].rp |> skipmissing |> collect
    s=  scatter(real.(samples).+randn(size(samples,1)).*0.1, 
        imag.(samples)+randn(size(samples,1)).*0.1, ms=1,
        label=""
    )
    X=[0+0im; mean(sp.rp|>skipmissing|>collect)]
    # Y = imag(X);
    # X = real(X); 
    plot!(real(X), imag(X), color=:black, label="$(nrow(sp)) spikes))")
    scatter!([real(X)[end]], [imag(X)[end]], color=:black, label="", ms=1)
    s
end
P = []
unit = first(unique(SPIKES.unit))
for unit in unique(SPIKES.unit)
    push!(P,plot_spike_phases(SPIKES[SPIKES.unit.==unit,:]))
end
plot(P[1:50]..., size=(1000,1000))
    
trp = combine(groupby(SPIKES, :unit, sort=true),
    :tp=>(x->skipmissing(x)|> collect |> length)=>:Nt,
    :rp=>(x->skipmissing(x)|> collect |> length)=>:Nr,
    :tp=>(x->skipmissing(x)|> collect |>mean)=>:tp,
    :rp=>(x->skipmissing(x)|>collect  |>mean)=>:rp,
                      renamecols=false)
trp = transform(trp, 
[:tp, :Nt] => ((tp, N)-> (N.*abs.(tp).^2 .- 1)./(N.-1)) => :tpc,
[:rp, :Nr] => ((rp, N)-> (N.*abs.(rp).^2 .- 1)./(N.-1)) => :rpc,
)
sort!(SPIKES, [:unit, :time])
DIutils.filtreg.register(trp, SPIKES, on="unit", transfer=["tpc", "rpc"])

histogram(angle.(disallowmissing(skipmissing(SPIKES.tp)|>collect)), 
    bins=100, normed=true, label="", 
    margins=1*Plots.mm,
    ylabel="Cells",
    title="Theta Phase",xlabel=L"Phase")

histogram(angle.(disallowmissing(skipmissing(SPIKES.rp)|>collect)), 
    bins=100, normed=true, label="", 
    margins=1*Plots.mm,
    ylabel="Cells",
    title="Ripple Phase",xlabel=L"Phase")

@assert !all(ismissing.(trp.rp))

# ================================================================
# Examining the theta and ripple phase modulation
# ================================================================

# Description: distribution of modulations
plot(
    (histogram(trp.tpc, bins=100, normed=true, label="", 
        margins=1*Plots.mm,
        ylabel="Cells",
        title="Theta Phase Modulation",xlabel=L"Modulation Index $\hat{\Gamma}$");
    vline!([0.01], label="",c=:black);),
    (histogram(trp.rpc, bins=100, normed=true, label="", 
        margins=1*Plots.mm,
        ylabel="Cells",
        xlabel=L"Modulation Index $\hat{\Gamma}$", title="Ripple Phase Modulation");
    vline!([0.02], label="",c=:black);),
)

# Description : how real and imag of ripple phase modulation relate to
#               corrected modulation
plot(
(scatter(real(trp.rp), trp.rpc, xlabel="real phase", ylabel="rip modulation");
    vline!([0], label="",c=:black);hline!([0], label="", c=:black)),
(scatter(imag(trp.rp), trp.rpc; xlabel="complex phase", ylabel="rip modulation");
    vline!([0], label="",c=:black);hline!([0], label="", c=:black)),
    title="Ripple Phase Modulation", label=""
)
# Description : how real and imag of theta phase modulation relate to
#              corrected modulation
plot(
(scatter(real(trp.tp), trp.tpc, xlabel="real phase", ylabel="theta modulation");
    vline!([0], label="",c=:black);hline!([0], label="", c=:black)),
(scatter(imag(trp.tp), trp.tpc; xlabel="complex phase", ylabel="theta modulation");
    vline!([0], label="",c=:black);hline!([0], label="", c=:black)),
    title="Theta Phase Modulation", label=""
)
# Description : how real and imag of theta and ripple phase modulation relate
# to each other
plot(
(scatter(real(trp.tp), real(trp.rp), xlabel="real theta phase",
        ylabel="real ripple phase");
    vline!([0], label="",c=:black);hline!([0], label="", c=:black)),
(scatter(imag(trp.tp), imag(trp.rp); xlabel="complex theta phase", 
        ylabel="complex ripple phase");
    vline!([0], label="",c=:black);hline!([0], label="", c=:black)),
    title="Phase Relationship", label=""
)

# Description : each cells Z mean arrow length modulation visualized
#               for theta
plot(
    (scatter(real(trp.tp), imag(trp.tp), projection=:polar, 
    xlabel="real theta phase", ylabel="complex theta phase", 
        title="Theta Phase", label="", markersize=3, alpha=0.5);
    vline!([0], c=:black, linestyle=:dash, label="");
    hline!([0], c=:black, linestyle=:dash, label="")
    ),
# Description : each cells Z mean arrow length modulation visualized
#               for ripple
    (scatter(real(trp.rp), imag(trp.rp), projection=:polar, 
    xlabel="real ripple phase", ylabel="complex ripple phase", 
    title="Ripple Phase", label="",
    markersize=3, alpha=0.5);
    vline!([0], c=:black, linestyle=:dash, label="");
    hline!([0], c=:black, linestyle=:dash, label="")
    )
)

# Create a 2x2 grid of subplots, each with a different camera angle
# First subplot
s1=scatter(real.(trp.rp), imag.(trp.rp), trp.rpc, camera=(0, 30), xlabel="real phase", ylabel="complex phase", zlabel="modulation")
# Second subplot
s2=scatter(real.(trp.rp), imag.(trp.rp), trp.rpc, camera=(45, 30), xlabel="real phase", ylabel="complex phase", zlabel="modulation")
# Third subplot
s3=scatter(real.(trp.rp), imag.(trp.rp), trp.rpc, camera=(90, 30), xlabel="real phase", ylabel="complex phase", zlabel="modulation")
# Fourth subplot
s4=scatter(real.(trp.rp), imag.(trp.rp), trp.rpc, camera=(135, 30), xlabel="real phase", ylabel="complex phase", zlabel="modulation")
# Show the plot
plot(s1, s2, s3, s4, layout=(2, 2), size=(800, 600), markersize=2)

#  . . . . . .FILTRATION . . . . . . . . . . . . . . . . . . . . . . . . . .
ca1units  = subset(CELLS, :area=>a->a.=="CA1", view=true).unit
pfcunits  = subset(CELLS, :area=>a->a.=="PFC", view=true).unit
SPIKESca1 = subset(SPIKES, :unit=>u->u.∈(ca1units,), view=true)
SPIKESpfc = subset(SPIKES, :unit=>u->u.∈(pfcunits,), view=true)
sp = SPIKESca1
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

# ------------------------------------------------------------------------
sptheta = subset(sp, :theta=>t->t  .!== missing .&& t.>0, 
                      :cycle => c->c.!==missing,
                     view=true)
sprip   = subset(sp, :ripple=>t->t .!== missing .&& t.>0,
                     :ripple => r-> r.!= NaN, view=true)
spthrip = subset(sp, :theta=>t->t .!== missing .&& t.>0, 
                     :ripple => r-> r.!= NaN, 
                     :ripple=>t->t .!== missing .&& t.>0, view=true)
# ------------------------------------------------------------------------

gr()

# ------------------------------------------------------------------------
# THETA AND RIPPLE PHASE DISTRIBUTIONS
# ------------------------------------------------------------------------
plot(
histogram(filter(t->t!=0, sptheta.theta_phase.*2pi), title="Theta phase distribution"),
histogram(filter(r->r!=0,sprip.ripple_phase*2pi),   title="Ripple phase distribution")
)
histogram2d(spthrip.ripple_phase, 
            spthrip.theta_phase, markersize=0.1, alpha=0.6, colorbar=true,
    xlabel="Ripple phase", ylabel="Theta phase", title="Theta-ripple phase relationship"
)
scatter!(spthrip.ripple_phase, spthrip.theta_phase, markersize=0.3333, alpha=0.6)

# ------------------------------------------------------------------------
# MEAN PHASE VERSUS TAU
# ------------------------------------------------------------------------
th = sort(dropmissing(combine(groupby(sptheta, :unit), :theta_phase=>mean,
    taus.=>first, renamecols=false), "all-1of20"),:theta_phase)
rip = sort(dropmissing(combine(groupby(sprip, :unit), :ripple_phase=>mean,
    taus.=>first, renamecols=false), "all-1of20"),:ripple_phase)
thrip = sort(dropmissing(combine(groupby(spthrip, :unit), :ripple_phase=>mean,
    :theta_phase=>mean, taus.=>first, renamecols=false), "all-1of20"),:ripple_phase)

# Add modulation scores
sort!(th, :unit)
DIutils.filtreg.register(trp, th, on="unit", transfer=["tpc"])
sort!(rip, :unit)
DIutils.filtreg.register(trp, rip, on="unit", transfer=["rpc"])

# Filter out units with low modulation
k=1
th  = filter(r->r.tpc>0.01*k, th)
rip = filter(r->r.rpc>0.01*k, rip)
thrip = filter(r->r.tpc>0.01*k && r.rpc>0.01*k, thrip)
theta_color  = :red
ripple_color = :skyblue

plot(
scatter(th.theta_phase,   th[!,"all-1of20"],  markersize=3, c=theta_color,
        xlabel="Theta Phase", ylabel="Tau", xlim=(0,1)),
scatter(rip.ripple_phase, rip[!,"all-1of20"], markersize=3, xlabel="Rip Phase",
        c=ripple_color,
        ylabel="Tau", xlim=(0,1)),
scatter(thrip.ripple_phase, thrip.theta_phase, thrip[!,"all-1of20"],
        markersize=3, xlabel="Ripple phase", ylabel="Theta phase",
        zlabel="Tau", xlim=(0,1)),
scatter(thrip.ripple_phase, thrip.theta_phase, markersize=3, alpha=0.5,
        xlim=(0,1), ylim=(0,1), xlabel="Ripple phase", ylabel="Theta phase",
    ),
title="Mean phase versus tau", size=(800,800)
    )

plot(
scatter(th.theta_phase,   th[!,"arena-1of20"],  markersize=3, c=theta_color,
        xlabel="Theta Phase", ylabel="Tau", xlim=(0,1)),
scatter(rip.ripple_phase, rip[!,"arena-1of20"], markersize=3, xlabel="Rip Phase",
        c=ripple_color,
        ylabel="Tau", xlim=(0,1)),
scatter(thrip.ripple_phase, thrip.theta_phase, thrip[!,"arena-1of20"],
        markersize=3, xlabel="Ripple phase", ylabel="Theta phase",
        zlabel="Tau", xlim=(0,1)),
scatter(thrip.ripple_phase, thrip.theta_phase, markersize=3, alpha=0.5,
        xlim=(0,1), ylim=(0,1), xlabel="Ripple phase", ylabel="Theta phase",
    ),
title="ARENA"
    )

plot(
scatter(th.theta_phase,   th[!,"home-1of20"],  markersize=3, c=theta_color,
        xlabel="Theta Phase", ylabel="Tau", xlim=(0,1)),
scatter(rip.ripple_phase, rip[!,"home-1of20"], markersize=3, xlabel="Rip Phase",
        c=:ripple_color,
        ylabel="Tau", xlim=(0,1)),
scatter(thrip.ripple_phase, thrip.theta_phase, thrip[!,"home-1of20"],
        markersize=3, xlabel="Ripple phase", ylabel="Theta phase",
        zlabel="Tau", xlim=(0,1)),
scatter(thrip.ripple_phase, thrip.theta_phase, markersize=3, alpha=0.5,
        xlim=(0,1), ylim=(0,1), xlabel="Ripple phase", ylabel="Theta phase",
    ),
title="home"
    )


#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
# Figure out correlated splits who could be examined together
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

nr = replace.(names(responses), "-1of20"=>"", "_" => " ")

# Raw heatmap
responses = th[!, Not([:unit, :theta_phase])]
Q = cor(Matrix(responses))
xticks = (1:size(Q,1), nr)
heatmap(Q, title="Correlation of TAU responses", clim=(-1,1), c=:vik,
    xticks=xticks, yticks=xticks, xrotation=45, yrotation=45,
    xtickfont=font(6), ytickfont=font(6), size=(800,800)
)

# Let's do a dendrogram in seaborn
using PyCall
il  = pyimport("importlib")
sns = pyimport("seaborn")
plt = pyimport("matplotlib.pyplot")
pd  = pyimport("pandas")
il.reload(sns); il.reload(plt)
plt.ion()
d = pd.DataFrame(Q, index=nr, columns=nr)
g = sns.clustermap(d, center=0, cmap="vlag", vmin=-1, vmax=1, 
    xticklabels=true, yticklabels=true, 
)
g.fig.show()
plt.pause(1)

g=nothing
plt.close("all")

# ------------------------------------------------------------------------
# TAU VERSUS THETA PHASE
# - home : check home
# - arena : check arena
# ------------------------------------------------------------------------
function phasetaudists(th; n="theta")
    phase = Symbol(n*"_phase")
    sort!(th, phase)
    responses = th[!, Not([:unit, phase])]
    nr = replace.(names(responses), "-1of20"=>"", "_" => " ")
    homes   = occursin.("home", nr)
    arenas  = occursin.("arena", nr)
    nontask = occursin.("nontask", names(responses))
    error = occursin.("error", names(responses))
    homes .&= .!error
    arenas .&= .!error
    homes   = names(responses)[homes]
    arenas  = names(responses)[arenas]
    nontask = names(responses)[nontask]
    error   = names(responses)[error]
    P = []
    for (filt, name) in zip([homes, arenas, nontask, error],
                            ["home", "arena", "nontask", "error"])
        TH = th[!, filt]
        responses = Matrix(TH)
        phases = repeat(th[!, phase], 1, size(responses,2))
        df = DataFrame(phases=phases[:], responses=responses[:])
        dfgt = dfgp = nothing
        begin
            df.phasebin = DIutils.binning.digitize(df.phases, 5)
            df.taubin  = DIutils.binning.digitize(disallowmissing(df.responses), 5)
            dfgt = combine(
                         groupby(df, :taubin), 
                         :phases=>mean,
                         :responses=>mean, renamecols=false)
            dfgp = combine(groupby(df, :phasebin), 
                         :phases=>mean,
                         :responses=>mean, renamecols=false)
        end
        begin
            df.phasebin = DIutils.binning.digitize(df.phases, 5)
            df.taubin  = DIutils.binning.digitize(disallowmissing(df.responses), 5)
            dfgt = combine(groupby(df, :taubin), 
                         :phases=>mean,
                         :responses=>mean, renamecols=false)
            dfgp = combine(groupby(df, :phasebin), 
                         :phases=>mean,
                         :responses=>mean, renamecols=false)
        end
        histogram2d(phases, responses, ylim=(-2,2))
        scatter!(phases, responses, markersize=1.25, alpha=1, colorbar=true,
            xlabel="$n phase", ylabel=L"$\tau$ timeshift", 
            title="$(uppercase(name))\n" * L" phase vs. $\tau$" * "\n", 
            label="", c=:darkgray, margin=5Plots.mm)
        plot!(dfgt.phases, dfgt.responses, linewidth=3, color=:black, 
            label="Mean", lw=3, marker=:circle, markersize=5, linestyle=:dash)
        plot!(dfgp.phases, dfgp.responses, linewidth=3, color=:black, 
            label="Mean", lw=3, marker=:circle, markersize=5, linestyle=:dash)
        vline!([mean(df.phases)], color=:black,     linestyle=:dash, label="Mean phase")
        hline!([mean(df.responses)], color=:black,  linestyle=:dash, label="Mean tau")
        push!(P, current())
    end
    plot(P..., size=(800,800), layout=(2,2), legend=:bottomright, 
        legendfontsize=3, 
    )
end
phasetaudists(th)
phasetaudists(rip;n="ripple")
