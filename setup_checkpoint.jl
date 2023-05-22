@time using GoalFetchAnalysis
import DI, DIutils
using DataFrames, Statistics
using Plots
using Infiltrator
animal, day = "super_clean", 0
@time SPIKES, BEH, CELLS, RIPPLES = DI.load(animal, day, 
    data_source=["spikes","behavior", "cells", "ripples"])
beh2 = DI.load_behavior("super",0)
SPIKES.cycle = Vector{Union{Float32, Missing}}(missing, size(SPIKES,1)) 
DI.annotate_pyrlayer!(CELLS)

Plot.setparentfolder("isolated")
# subset(combine(groupby(cells, [:animal,:tetrode,:pyrlayer]),
#     :meanrate=>mean,
#     :meanrate=>length=>:n_cells,
#     :propbursts=>mean,
# ), :pyrlayer=>s->s.!=false)
# Time conversion factors
# time_factors = DI.superanimal_timeconversion("super",0)
time_factors  = DI.behavior_0time_before_process("super")
animals, days = unique(CELLS.animal), unique(CELLS.day)
animal, day   = zip(animals, days)|>first
DIutils.pushover("Finished loading ALL ANIMAL data")
animal = animals[1]
load_from_checkpoint = true

        l_pyr = DI.load_lfp(animal, day; append="pyr")
        l_pyr = transform(l_pyr, :time=> t-> t .- time_factors[animal], renamecols=false)
        # transform!(l_pyr, :time=> t-> t .+ time_factors[animal], renamecols=false)
        cycles       = DI.load_cycles(animal, day, "pyr")
        cycles       = transform(cycles,
            :start=> t-> t .- time_factors[animal],
            :stop=> t-> t .- time_factors[animal],
            renamecols=false
        )
        spikes = transform(DI.load_spikes(animal, day, "pyr_cycles_isolated"),
            :time=> t-> t .- time_factors[animal],
            renamecols=false
        )
        # spikes = subset(SPIKES, :animal=>a->a.==animal, :day=>d->d.==day)
        # Print out extrema(time) for each of these just to be sure
        
        begin
            println("spikes.time: ", extrema(spikes.time))
            println("lfp.time: ", extrema(l_pyr.time))
            println("cycles.time: ", extrema(cycles.start))
        end
        DIutils.pushover("Loaded checkpoint for $animal")

convert_to_f32 = [:broadraw  :phase :ripple  :rippleamp  :ripplephase]
for col in convert_to_f32
    if hasproperty(l_pyr, col)
        l_pyr[!,col] = convert(Array{Float32,1}, l_pyr[!,col])
    end
end
