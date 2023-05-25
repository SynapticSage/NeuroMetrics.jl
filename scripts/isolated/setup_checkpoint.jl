@time using GoalFetchAnalysis
import DI, DIutils
using DataFrames, Statistics
using Plots
using Infiltrator
function checkranges()
    if isdefined(Main,:spikes)
        println("Spikes.time extrema: ", extrema(spikes.time))
        @assert minimum(spikes.time) <= 0
    end
    if isdefined(Main,:l_pyr)
        println("l_pyr.time extrema: ", extrema(l_pyr.time))
        @assert minimum(l_pyr.time) <= 0
    end
    if isdefined(Main,:cycles)
        println("cycles.time extrema: ", extrema(cycles.start))
        @assert minimum(cycles.start) <= 0
    end
    if isdefined(Main,:lfp) && lfp !== nothing
        println("lfp.time extrema: ", extrema(lfp.time))
        @assert minimum(lfp.time) <= 0
    end
    if isdefined(Main,:ripple) && ripple !== nothing
        println("ripple.time extrema: ", extrema(ripple.time))
        @assert minimum(ripple.time) <= 0
    end
end
checkranges()

function checkanimalranges()
    println("BEH:\n",
        combine(groupby(BEH, [:animal, :day]),    :time=>minimum),
        "\nSPIKES:\n",
            combine(groupby(SPIKES, [:animal, :day]), :time=>minimum))
end

# load super animal data
animal, day = "super_clean", 0
@time SPIKES, BEH, CELLS, RIPPLES = DI.load(animal, day, 
    data_source=["spikes","behavior", "cells", "ripples"])
beh2 = DI.load_behavior("super",0)
SPIKES.cycle = Vector{Union{Float32, Missing}}(missing, size(SPIKES,1)) 
DI.annotate_pyrlayer!(CELLS)
Plot.setparentfolder("isolated")
lfp_convert    = DI.get_0time_pre_superanimal()
super_convert  = DI.convert_super_to_time0()
load_pyr=nothing
if !isdefined(Main, :opt)
    animals, days = unique(CELLS.animal), unique(CELLS.day)
    animal, day   = collect(zip(animals, days))[2]
    load_pyr = true
else
    animal, day = opt["animal"], opt["day"]
    load_pyr = opt["load_pyr"]
end
opt["animal"], opt["day"], opt["load_pyr"] = animal, day, load_pyr
if load_pyr
    println("Loading PYR")
    l_pyr = DI.load_lfp(animal, day; append="pyr")
    l_pyr = transform(l_pyr, :time=> t-> t .- lfp_convert[animal], renamecols=false)
    # transform!(l_pyr, :time=> t-> t .+ time_factors[animal], renamecols=false)
    cycles       = DI.load_cycles(animal, day, "pyr")
    cycles       = transform(cycles,
        :start=> t-> t .- lfp_convert[animal],
        :stop=> t-> t .- lfp_convert[animal],
        renamecols=false
    )
    spikes = transform(DI.load_spikes(animal, day, "pyr_cycles_isolated"),
        :time=> t-> t .- lfp_convert[animal],
        renamecols=false
    )
    spikes = subset(SPIKES, :animal=>a->a.==animal, :day=>d->d.==day)
    # Print out extrema(time) for each of these just to be sure
    begin
        checkranges()
    end
    convert_to_f32 = [:broadraw  :phase :ripple  :rippleamp  :ripplephase]
    for col in convert_to_f32
        if hasproperty(l_pyr, col)
            l_pyr[!,col] = convert(Array{Float32,1}, l_pyr[!,col])
        end
    end
else
    println("lfp_convert: ", lfp_convert)
    println("super_convert: ", super_convert)
    checkanimalranges()
    beh = subset(BEH, :animal=>a->a.==animal, :day=>d->d.==day)
    for animal in animals
        println("Correcting $animal by ", super_convert[animal])
        BEH[BEH.animal.==animal, :time] .-= super_convert[animal]
        SPIKES[SPIKES.animal.==animal, :time] .-= super_convert[animal]
    end
    checkanimalranges()
    spikes = subset(SPIKES, :animal=>a->a.==animal, :day=>d->d.==day)
end

checkranges()
DIutils.pushover("Finished loading $animal")

