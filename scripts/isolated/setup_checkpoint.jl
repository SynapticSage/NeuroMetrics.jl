@time using GoalFetchAnalysis
import DI, DIutils
using DataFrames, Statistics
using Plots
using DataFrames
using Infiltrator
using Peaks

function checkranges()
    if isdefined(Main,:spikes) && spikes !== nothing
        println("Spikes.time extrema: ", extrema(spikes.time))
        # @assert minimum(spikes.time) <= 0
    else
        @warn "spikes is not defined"
    end
    if isdefined(Main,:cycles) && cycles !== nothing
        println("cycles.time extrema: ", extrema(cycles.start))
        # @assert minimum(cycles.start) <= 0
    else
        @warn "cycles is not defined"
    end
    if isdefined(Main,:lfp) && lfp !== nothing
        println("lfp.time extrema: ", extrema(lfp.time))
        # @assert minimum(lfp.time) <= 0
    else
        @warn "lfp is not defined"
    end
    if isdefined(Main,:ripples) && ripples !== nothing
        println("ripple.time extrema: ", extrema(ripples.time))
        # @assert maximum(ripples.time) > 0
    else
        @warn "ripples is not defined"
    end
end

function checkanimalranges()
    println(
        "BEH:\n",
        combine(groupby(BEH, [:animal, :day]),    :time=>minimum),
        "\nSPIKES:\n",
            combine(groupby(SPIKES, [:animal, :day]), :time=>minimum),
        "\nRIPPLES:\n",
        isdefined(Main, :RIPPLES) ?
        combine(groupby(RIPPLES, [:animal, :day]), :time=>minimum) :
                            nothing
    )
end

# -----------------------
# load super animal data
# -----------------------
animal, day = "super_clean", 0
@time SPIKES, BEH, CELLS = DI.load(animal, day, 
    data_source=["spikes","behavior", "cells"])
beh2 = DI.load_behavior("super",0)
SPIKES.cycle = Vector{Union{Float32, Missing}}(missing, size(SPIKES,1)) 
DI.annotate_pyrlayer!(CELLS)
Plot.setparentfolder("isolated")
lfp_convert    = DI.get_0time_pre_superanimal()
super_convert  = DI.convert_super_to_time0()
load_pyr=nothing
animals, days = unique(CELLS.animal), unique(CELLS.day)
if !isdefined(Main, :opt)
    animal, day   = collect(zip(animals, days))[2]
    load_pyr = true
else
    animal, day = opt["animal"], opt["day"]
    load_pyr    = opt["load_pyr"]
end


# -----------------------
# Load specific animal data
# -----------------------
if !isdefined(Main, :opt)
    # opt = Dict(
    #     "animal"=>animal,
    #     "day"=>day,
    #     "load_pyr"=>load_pyr
    # )
    opt = Dict(
        "animal"=>"RY16",
        "day"=>day,
        "load_pyr"=>false,
    )
end

function get_animal_pyr(animal, day)
    println("Loading PYR data for $animal")
    lfp = DI.load_lfp(animal, day; append="$tetrode_set")
    lfp = transform(lfp, :time=> t-> t .- lfp_convert[animal],
        renamecols=false)
    # transform!(lfp, :time=> t-> t .+ time_factors[animal], renamecols=false)
    cycles       = DI.load_cycles(animal, day, "pyr")
    cycles       = transform(cycles,
        :start=> t-> t .- lfp_convert[animal],
        :stop=> t-> t .- lfp_convert[animal],
        renamecols=false
    )
    global spikes = transform(DI.load_spikes(animal, day, "pyr_cycles_isolated"),
        :time=> t-> t .- lfp_convert[animal],
        renamecols=false
    )
    spikes = subset(SPIKES, :animal=>a->a.==animal, :day=>d->d.==day)
    # Print out extrema(time) for each of these just to be sure
    convert_to_f32 = [:broadraw  :phase :ripple  :rippleamp  :ripplephase]
    for col in convert_to_f32
        if hasproperty(lfp, col)
            lfp[!,col] = convert(Array{Float32,1}, lfp[!,col])
        end
    end
    global ripples = DI.load_ripples(animal, day)
    ripples = transform(ripples, 
        :time=> t-> t .- lfp_convert[animal],
        :start=> t-> t .- lfp_convert[animal],
        :stop=> t-> t .- lfp_convert[animal],
        renamecols=false)
    return lfp, cycles, spikes, ripples
end

# MULTIPLE ANIMALS
if load_pyr && !occursin("super",animal)
    global lfp, cycles, spikes, ripples = get_animal_pyr(animal, day)
    checkranges()
# ONE ANIMAL
elseif load_pyr
    OUT = []
    animals, days = DI.animaldays()
    for (animal,day) in zip(animals, days)
        lfp, cycles, spikes, ripples = get_animal_pyr(animal, day)
        push!(OUT, (lfp, cycles, spikes, ripples))
    end
    global lfp = vcat([x[1] for x in OUT]...)
    global cycles = vcat([x[2] for x in OUT]...)
    global spikes = vcat([x[3] for x in OUT]...)
    global ripples = vcat([x[4] for x in OUT]...)
# ONE ANIMAL from scratch
else
    println("lfp_convert: ", lfp_convert)
    println("super_convert: ", super_convert)
    checkanimalranges()
    beh = subset(BEH, :animal=>a->a.==animal, :day=>d->d.==day)
    for animal in animals
        println("Correcting $animal by ", super_convert[animal])
        BEH[BEH.animal         .== animal, :time] .-= super_convert[animal]
        SPIKES[SPIKES.animal   .== animal, :time] .-= super_convert[animal]
    end
    checkanimalranges()
    global spikes = subset(SPIKES, :animal=>a->a.==animal, :day=>d->d.==day, 
                     view=true)
    global ripples = DI.load_ripples(animal, day)
    global ripples = transform(ripples, 
        :time=> t-> t .- lfp_convert[animal],
        :start=> t-> t .- lfp_convert[animal],
        :stop=> t-> t .- lfp_convert[animal],
        renamecols=false)
end

checkranges()
DIutils.pushover("Finished loading $animal")

