@time using GoalFetchAnalysis
import GoalFetchAnalysis: Munge
import GoalFetchAnalysis: Plot
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
beh2 = DI.load_behavior("super_clean",0)
SPIKES.cycle = Vector{Union{Float32, Missing}}(missing, size(SPIKES,1)) 
DI.annotate_pyrlayer!(CELLS)
Plot.setparentfolder("isolated")
lfp_convert    = DI.get_0time_pre_superanimal()
super_convert  = DI.convert_super_to_time0()
checkpoint=nothing
animals, days = unique(CELLS.animal), unique(CELLS.day)

if !isdefined(Main, :opt)
    animal, day = collect(zip(animals, days))[2]
    checkpoint  = true
else
    animal, day = opt["animal"], opt["day"]
    checkpoint  = opt["checkpoint"]
end


# -----------------------
# Load specific animal data
# -----------------------
if !isdefined(Main, :opt)
    # opt = Dict(
    #     "animal"=>animal,
    #     "day"=>day,
    #     "checkpoint"=>checkpoint
    # )
    opt = Dict(
        "animal"=>"RY16",
        "day"=>day,
        "checkpoint"=>false,
    )
end
cells, spikes, beh, ripples, cycles, lfp = 
    Munge.isolated.import_dataframes(opt["animal"], opt["day"]; checkpoint=opt["checkpoint"])

checkranges()
DIutils.pushover("Finished loading $animal")

