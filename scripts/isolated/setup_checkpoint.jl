@time using GoalFetchAnalysis
import DI, DIutils
using DataFrames, Statistics
using Plots
using Infiltrator
# load super animal data
animal, day = "super_clean", 0
@time SPIKES, BEH, CELLS, RIPPLES = DI.load(animal, day, 
    data_source=["spikes","behavior", "cells", "ripples"])
beh2 = DI.load_behavior("super",0)
SPIKES.cycle = Vector{Union{Float32, Missing}}(missing, size(SPIKES,1)) 
DI.annotate_pyrlayer!(CELLS)
Plot.setparentfolder("isolated")
time_factors  = DI.behavior_0time_before_process("super")
load_pyr=nothing
if !isdefined(Main, :opt)
    animals, days = unique(CELLS.animal), unique(CELLS.day)
    animal, day   = collect(zip(animals, days))[2]
    load_pyr = true
else
    animal, day = opt["animal"], opt["day"]
    load_pyr = opt["process_lfp"]
end
# if load_pyr
#     println("Loading PYR")
#     l_pyr = DI.load_lfp(animal, day; append="pyr")
#     l_pyr = transform(l_pyr, :time=> t-> t .- time_factors[animal], renamecols=false)
#     # transform!(l_pyr, :time=> t-> t .+ time_factors[animal], renamecols=false)
#     cycles       = DI.load_cycles(animal, day, "pyr")
#     cycles       = transform(cycles,
#         :start=> t-> t .- time_factors[animal],
#         :stop=> t-> t .- time_factors[animal],
#         renamecols=false
#     )
#     spikes = transform(DI.load_spikes(animal, day, "pyr_cycles_isolated"),
#         :time=> t-> t .- time_factors[animal],
#         renamecols=false
#     )
#     spikes = subset(SPIKES, :animal=>a->a.==animal, :day=>d->d.==day)
#     # Print out extrema(time) for each of these just to be sure
#     begin
#         println("spikes.time: ", extrema(spikes.time))
#         println("lfp.time: ", extrema(l_pyr.time))
#         println("cycles.time: ", extrema(cycles.start))
#     end
#     convert_to_f32 = [:broadraw  :phase :ripple  :rippleamp  :ripplephase]
#     for col in convert_to_f32
#         if hasproperty(l_pyr, col)
#             l_pyr[!,col] = convert(Array{Float32,1}, l_pyr[!,col])
#         end
#     end
# end
DIutils.pushover("Finished loading $animal")

function checkranges()
begin
    if isdefined(Main,:spikes)
        println("Spikes.time extrema: ", extrema(spikes.time))
    end
    if isdefined(Main,:l_pyr)
        println("l_pyr.time extrema: ", extrema(l_pyr.time))
    end
    if isdefined(Main,:cycles)
        println("cycles.time extrema: ", extrema(cycles.start))
    end
    if isdefined(Main,:lfp)
        println("lfp.time extrema: ", extrema(lfp.time))
    end
end
end
checkranges()
