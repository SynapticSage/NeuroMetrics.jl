### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook
# outside of Pluto, the following 'mock version' of @bind gives bound variables
# a default value (instead of an error).
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
	import Plot
	Field, Load, Timeshift, Utils, Load, Table, Plot
end

# ╔═╡ 81dd0cb9-5400-4502-8a27-0530d730c616
begin
	import PlutoUI
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
	using DataStructures

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

# ╔═╡ 94816c70-20ab-4e80-93fb-ed1ce0b3405a
PlutoUI.TableOfContents(title="Comparing field shifts")

# ╔═╡ 6a1d85d7-f98e-4af5-aa8d-0f9696b8bad8
md"""
# Preamble

## Imports
"""

# ╔═╡ 4d85a1f8-fe28-45bd-970c-cab78838034f
md"""
## Load data

Load up shifting data and subset it parametrically.
"""

# ╔═╡ 2ab79319-1f6f-431f-93ab-dcb5c381afd8
I, cells = with_logger(NullLogger()) do
cells = Load.load_cells("RY16", 36, "*");
_, I = Load.register(cells, Timeshift.load_mains(dataframe=true), on="unit", transfer="meanrate");
	I, cells
end;

# ╔═╡ 68f79249-c6f2-4f54-b3b6-e0559f9bb7b5
md"""
### Parameters
"""

# ╔═╡ 49896b0c-19d3-4c51-8e54-2f046e6d8382
md"""
We have some choices to make. This worksheet will compare any *two* checked datasets.  In the future, maybe N. But for now, 2 is the limit.
"""

# ╔═╡ e8221624-f6fd-11ec-0151-cf9da8e7d639
begin
    dataset₁ = @bind datacut₁ PlutoUI.MultiSelect(String.(unique(I.datacut)), default=["all"])
	 dataset₂ = @bind datacut₂ PlutoUI.MultiSelect(String.(unique(I.datacut)), default=["all"])
	dataset₁, dataset₂
end

# ╔═╡ 948387de-df9e-4dde-b616-4fa54e13738d
datacut = [datacut₁..., datacut₂...]

# ╔═╡ b38a32ff-9b95-416c-be82-3d7371d51932
@bind marginal PlutoUI.Radio(String.(unique(I.marginal)), default="x-y")

# ╔═╡ b741339a-8a3b-416d-9a59-33624595f56f
@bind area PlutoUI.Radio(["CA1-PFC", String.(unique(I.area))...], default="CA1-PFC")

# ╔═╡ 4c670dcd-25ff-43d1-b78f-2b60fc5878d5
comparison_name = "Δ $(datacut[1])-$(datacut[2])"

# ╔═╡ 48c22294-fc9d-44ff-b2ea-b2c59e795da0
datacut

# ╔═╡ 61f417ed-a4c9-42a3-b858-17dc170dfea0
Is = begin
	Is = OrderedDict{Symbol,DataFrame}()
	for cut in datacut
		tmp = transform(DataFramesMeta.@subset(I, :datacut .== Symbol(cut), :marginal .== marginal), :shift => (x->x*60) => :shift)
		if area != "CA1-PFC"
			@subset!(tmp, :area .== area)
		end
		Is[Symbol(cut)] = tmp
	end
	
	TT = Table.group.multitable_groupby([:unit, :shift], values(Is)...);
	dfDiff = DataFrame()
	for (a,b) in eachrow(TT)
		c = copy(a)
		if ismissing(a) || ismissing(b)
			continue
		end
		
		c.datacut .= Symbol("compare_" * 
							String(a.datacut[1]) * "_" * 
							String(b.datacut[1]))
		c.value = a.value - b.value
		append!(dfDiff, c)
	end
	dfDiff=Timeshift.operation.add_best_shift(dfDiff);
	dfDiff=Timeshift.operation.add_centroid_shift(dfDiff);
	dfDiff
end;

# ╔═╡ a7727c6b-b284-4862-876d-591051ce5a2d
with_logger(NullLogger()) do
h_diff = heatmap(Float32.(replace(
	Matrix(
		sort(unstack(Timeshift.norm_information(dfDiff), [:unit, :bestshift], :shift, :value)[:,Not(:unit)], :bestshift)
	), missing=>NaN)
), clim=(0.25,1))
end

# ╔═╡ ee2ce4a3-f36f-4814-9853-1dddf1984b56
datasets

# ╔═╡ 11704bde-db87-4469-97cf-e2a4fc5c4c76
begin
	cs_diff = :oxy
	
	p1 = @df Timeshift.norm_information(dfDiff,minmax=[-1,1]) density(:value, xlabel=comparison_name,label="Overall difference")
	hline!([0.5], linestyle=:dash,label="")
	vline!([0], linestyle=:dash, c=:gray,label="")
	
	p2 = @df Timeshift.norm_information(dfDiff,minmax=[-1,1]) density(:value, group=:shift, legend=:none, palette=palette(cs_diff,41), xlim=(-1,1), title=String(dfDiff.datacut[1]) * " @ shifts", xlabel=comparison_name)
	hline!([0.5], linestyle=:dash, c=:gray)
	vline!([0], linestyle=:dash, c=:gray)

	cbar=heatmap([xlims(p2)...], [ylims(p2)...], repeat([NaN], 2,2), clim=(-4,4), colorbar=true,c=cs_diff, cbar_title="Shifts", yticks=[], xticks=[])
	layout = @layout [Plots.grid(1,2) b{0.075w}]
	plot(p1,p2,cbar; link=:y, layout)
end

# ╔═╡ Cell order:
# ╟─94816c70-20ab-4e80-93fb-ed1ce0b3405a
# ╟─6a1d85d7-f98e-4af5-aa8d-0f9696b8bad8
# ╟─ef3ed745-3204-494b-bd4a-f0c849255f2e
# ╟─f9d33aee-5392-41a6-bf03-fd9244adfa0f
# ╟─81dd0cb9-5400-4502-8a27-0530d730c616
# ╟─4d85a1f8-fe28-45bd-970c-cab78838034f
# ╟─2ab79319-1f6f-431f-93ab-dcb5c381afd8
# ╟─68f79249-c6f2-4f54-b3b6-e0559f9bb7b5
# ╟─49896b0c-19d3-4c51-8e54-2f046e6d8382
# ╟─e8221624-f6fd-11ec-0151-cf9da8e7d639
# ╟─948387de-df9e-4dde-b616-4fa54e13738d
# ╟─b38a32ff-9b95-416c-be82-3d7371d51932
# ╟─b741339a-8a3b-416d-9a59-33624595f56f
# ╟─4c670dcd-25ff-43d1-b78f-2b60fc5878d5
# ╟─48c22294-fc9d-44ff-b2ea-b2c59e795da0
# ╠═61f417ed-a4c9-42a3-b858-17dc170dfea0
# ╠═a7727c6b-b284-4862-876d-591051ce5a2d
# ╟─ee2ce4a3-f36f-4814-9853-1dddf1984b56
# ╠═11704bde-db87-4469-97cf-e2a4fc5c4c76
