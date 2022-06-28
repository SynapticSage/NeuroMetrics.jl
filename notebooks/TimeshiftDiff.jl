### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ ef3ed745-3204-494b-bd4a-f0c849255f2e
begin
	using Revise
	import DrWatson
	DrWatson.quickactivate(expanduser("~/Projects/goal-code"))
	push!(LOAD_PATH, DrWatson.srcdir())
end

# ╔═╡ f9d33aee-5392-41a6-bf03-fd9244adfa0f
Field, Load, Timeshift, Utils, Load, Table = begin
	using GoalFetchAnalysis
	import Timeshift
	import Field
	import Utils
	import Load
	import Table
	Field, Load, Timeshift, Utils, Load, Table
end

# ╔═╡ 81dd0cb9-5400-4502-8a27-0530d730c616
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

# ╔═╡ 6a1d85d7-f98e-4af5-aa8d-0f9696b8bad8
md"""
# Preamble

## Imports
"""

# ╔═╡ 4d85a1f8-fe28-45bd-970c-cab78838034f
md"""
## Load data
"""

# ╔═╡ e8221624-f6fd-11ec-0151-cf9da8e7d639
@bind datacut PlutoUI.Radio(String.(unique(I.datacut)), default="all")

# ╔═╡ b38a32ff-9b95-416c-be82-3d7371d51932
@bind marginal PlutoUI.Radio(String.(unique(I.marginal)), default="x-y")

# ╔═╡ 68f79249-c6f2-4f54-b3b6-e0559f9bb7b5
md"""
## Parameters
"""

# ╔═╡ Cell order:
# ╟─6a1d85d7-f98e-4af5-aa8d-0f9696b8bad8
# ╟─ef3ed745-3204-494b-bd4a-f0c849255f2e
# ╠═f9d33aee-5392-41a6-bf03-fd9244adfa0f
# ╠═81dd0cb9-5400-4502-8a27-0530d730c616
# ╠═4d85a1f8-fe28-45bd-970c-cab78838034f
# ╟─e8221624-f6fd-11ec-0151-cf9da8e7d639
# ╟─b38a32ff-9b95-416c-be82-3d7371d51932
# ╠═68f79249-c6f2-4f54-b3b6-e0559f9bb7b5
