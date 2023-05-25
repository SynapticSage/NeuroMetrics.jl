### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 20069172-a26b-449c-b6a7-36ae9f027fc3
begin
	begin
		import Pkg
		Pkg.activate(expanduser("~/Projects/goal-code"))
	end
end

# ╔═╡ cb4250bd-bb41-4aef-8f81-9e26bbdc5a95
begin
	using GoalFetchAnalysis
	using DataFrames, DataFramesMeta, ProgressMeter, StatsBase, Statistics,
	StatsPlots, Plots, LaTeXStrings, DirectionalStatistics
	using DSP: angle
	Base.angle(x::Missing) = missing
	using LaTeXStrings
	import Random, Pluto
end

# ╔═╡ f66c7e02-f8d9-11ed-2902-999801e69fce
begin
#=
=====================================================
Script Name     : VS_lfp.jl
Autho=r          : Ryan Young
Created         : 05-03-2023
Last Modified   : 05-03-2023
===================================================

# Description
This script plots local field potential properties against timeshifted
optimal 𝛕 of each cell.
=#
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
using Infiltrator
function get_data(animal, day)
    # l_pyr = DI.load_lfp(animal, day; append="pyr")
    # l_pyr = transform(l_pyr, :time=> t-> t .- time_factors[animal], renamecols=false)
    # transform!(l_pyr, :time=> t-> t .+ time_factors[animal], renamecols=false)
    cycles = DI.load_cycles(animal, day, "pyr")
    cycles = transform(cycles,
        :start=> t-> t .- time_factors[animal],
        :stop=> t-> t .-  time_factors[animal],
        renamecols=false
    )
    # cycles = cycles[!,Not([:start_function,:stop_function])]
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
end

# ╔═╡ 62a0b782-e556-4569-9aa3-93f26f1cd0ad
md"""# Section 1"""

# ╔═╡ 992e99bb-5ed2-40ea-88b4-917d7b42fd3d
md"""# Section 2"""

# ╔═╡ aae97da6-eeea-4ab3-819c-6073afd019ba
md"""# Section 3"""

# ╔═╡ 409237f5-a07a-4fe0-a186-2218fc62655c


# ╔═╡ Cell order:
# ╠═20069172-a26b-449c-b6a7-36ae9f027fc3
# ╠═cb4250bd-bb41-4aef-8f81-9e26bbdc5a95
# ╠═f66c7e02-f8d9-11ed-2902-999801e69fce
# ╠═62a0b782-e556-4569-9aa3-93f26f1cd0ad
# ╠═992e99bb-5ed2-40ea-88b4-917d7b42fd3d
# ╠═aae97da6-eeea-4ab3-819c-6073afd019ba
# ╠═409237f5-a07a-4fe0-a186-2218fc62655c
