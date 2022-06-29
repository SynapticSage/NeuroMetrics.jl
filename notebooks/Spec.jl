### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ dfcfa875-d122-49f1-ab24-66c1937b3134
begin
	using Revise
	import DrWatson
	DrWatson.quickactivate(expanduser("~/Projects/goal-code"))
	push!(LOAD_PATH, DrWatson.srcdir())
end

# ╔═╡ a9b4b3d2-f318-11ec-210a-a70a7964ee72
Field, Load, Timeshift, Utils, Load, Table, Plot, Munge = begin
	using GoalFetchAnalysis
	import Timeshift
	import Field
	import Utils
	import Load
	import Table
	import Plot
	import Munge
	Field, Load, Timeshift, Utils, Load, Table, Plot, Munge
end

# ╔═╡ 4db1e3a6-cc16-404c-beb8-0b2e60e19d59
using FourierAnalysis

# ╔═╡ 450738b2-3d49-4e45-9a4d-ffa1721f833a
begin
	using StatsBase
	using Term
	using REPL.TerminalMenus
	using Plots
	import Arrow
	using DataFramesMeta
	using DataFrames
	using ColorSchemes
	using StatsBase
	using Infiltrator
	using Dates
	using StatsPlots
	using Logging

	if Utils.in_range(hour(now()), [0,5]) ||
	   Utils.in_range(hour(now()), [20, 24])
		Plots.theme(:dark)
		theme="dark"
	else
		Plots.theme(:bright)
		theme="bright"
	end
	
	"Importing packages, theme=$theme"
end

# ╔═╡ 435f9680-1520-468e-b97c-2ea4fb2c1ff4
using PlutoUI

# ╔═╡ 8d41c178-16ee-4881-b55c-fb80f274d7bc
PlutoUI.TableOfContents(title="Non-adaptive field shifting")

# ╔═╡ 42ea762b-12ed-4eb8-ade0-3bffff593690
md"""
# Package Imports 📦
"""

# ╔═╡ bf904f0e-7387-4b73-890b-fbb5fc6137de
md"""
# Loading data  💾
First, we're going to loadup a key dataframe of interest: **cells**
"""

# ╔═╡ a7d10b7e-6534-4017-a042-1a414a9e96bb
md"""
# Spectra 🌊
 
Obtain the spectra
"""


# ╔═╡ 56b09cc7-e8c2-439c-9b02-2fb5dab44f72
spectra(lfp.raw)

# ╔═╡ cadaf555-3b90-4c5b-846b-686ce4130497
begin
	@time spikes, beh, ripples, cells = Load.load("RY16", 36);
	lfp = Load.load_lfp("RY16", 36)
end;

# ╔═╡ e99882e8-7a93-4569-b7ce-c4b83df58510
lfp = Munge.lfp.get_tet(sort(lfp, [:epoch,:time]), 5)

# ╔═╡ Cell order:
# ╠═8d41c178-16ee-4881-b55c-fb80f274d7bc
# ╟─dfcfa875-d122-49f1-ab24-66c1937b3134
# ╟─42ea762b-12ed-4eb8-ade0-3bffff593690
# ╟─a9b4b3d2-f318-11ec-210a-a70a7964ee72
# ╠═4db1e3a6-cc16-404c-beb8-0b2e60e19d59
# ╟─450738b2-3d49-4e45-9a4d-ffa1721f833a
# ╟─bf904f0e-7387-4b73-890b-fbb5fc6137de
# ╟─cadaf555-3b90-4c5b-846b-686ce4130497
# ╠═e99882e8-7a93-4569-b7ce-c4b83df58510
# ╠═a7d10b7e-6534-4017-a042-1a414a9e96bb
# ╠═56b09cc7-e8c2-439c-9b02-2fb5dab44f72
# ╟─435f9680-1520-468e-b97c-2ea4fb2c1ff4
