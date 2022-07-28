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

# ╔═╡ 44dde9e4-f9ca-11ec-1348-d968780f671c
begin
	  using DrWatson
	  quickactivate(expanduser("~/Projects/goal-code"))
	  using GoalFetchAnalysis
	  using Plots
	  using Revise
	  using DataFrames
	  props = ["x","y"]
	  using NaNStatistics
	  adaptive = Field.adaptive
	  import Utils
end

# ╔═╡ d33ca749-25af-4b1e-b78c-00ebf54f1327
using PlutoUI

# ╔═╡ 0be7ba01-a316-41ce-8df3-a5ae028c74e7
PlutoUI.TableOfContents(title="Test adaptive")

# ╔═╡ 37d7f4fd-80a7-47d0-8912-7f002620109f
md"""
Preamble ... loading and importing
"""

# ╔═╡ 31082fe7-ed61-4d37-a025-77420da3f24a
begin
	
	@time spikes, beh, ripples, cells = Load.load("RY16", 36);
	@time beh, spikes  = Load.register(beh, spikes; on="time", transfer=props)
	
end;

# ╔═╡ d51ce0f2-03bf-4c88-9302-9fd4bc8621eb
width_select = @bind width Slider(1f0:0.25f0:8f0, show_value=true, default=2.5f0)

# ╔═╡ efcdc2f1-5e26-4534-953e-defae4bd8603
md"""
# Grid
Let's try making a grid object
"""

# ╔═╡ 6bdcf863-9946-4ca3-ab02-fa6aebe4b91d
G = @time adaptive.get_grid(beh, props; widths=width)

# ╔═╡ 03586347-83ee-429d-ab29-505754c66734
plot(G)

# ╔═╡ 7c03a7aa-3181-4b3d-9dc6-0eb8e8689023
md"""
# Occupancy
And an occupancy object
"""

# ╔═╡ 890fe951-19bb-4c9a-a905-d798bb36c57e
O = @time adaptive.get_occupancy(beh, G)

# ╔═╡ 592d79b4-edf6-4a0c-af73-1d2805d6410e
plot(O, clim=(0,0.001))

# ╔═╡ 2c794dc4-4920-4274-ab2b-6bb60251112b
md"""
# Single adaptive field
grab all spikes as single multiunit field
"""

# ╔═╡ bef016cd-26d1-4de8-a970-182fe2b92e88
# Test field abilities
@time multiunit = @time adaptive.get_adaptivefield(spikes, G, O)

# ╔═╡ 410128bf-2332-4aa2-91cb-441e0235cc4a
plot(multiunit)

# ╔═╡ c23ee21a-b3d4-42e5-a4f1-b703008eda1a
md"""
# Field dict of all cells
"""

# ╔═╡ b88c0ec1-b150-49be-828f-6c32bb770c48
@time units = adaptive.ulanovsky(spikes, G, O; splitby=[:unit], widths=width)

# ╔═╡ 38cff24f-bbc1-42fd-98ae-385323c2480e
width_select

# ╔═╡ 7c6cfeb1-2c78-4480-852b-aa06cc818f76
unit_select = @bind unit PlutoUI.Slider(sort(unique(spikes.unit)), show_value=true)

# ╔═╡ 44abcbd4-5f71-4924-b77d-9680cc96044f
plot(units[(;unit=unit)], aspect_ratio=1)

# ╔═╡ 4d814c3e-97e1-491a-b1d8-c7ca9c628afd
μ_firing = begin
		Q = units[(;unit=unit)]
		nansum(reshape(Q.occ.prob, size(Q.occ.count)) .* Q.rate)
	end

# ╔═╡ f9378d49-2f86-4088-bc6d-3b5b227b7c66
md"""
## `to_dataframe`
"""

# ╔═╡ 34b5441c-add2-4272-b384-67994daf7745
md"""
# Timeshift
Getting everything working with Timeshift.jl
"""

# ╔═╡ 9e5fde4a-f312-419a-b7d5-5ef7a08c305e
import Timeshift

# ╔═╡ 39737bd9-f38a-408d-a0c0-99b9e2bd0045
shifted = Timeshift.shifted_fields(beh, spikes, -2:0.2:2, props, widths=3f0)

# ╔═╡ 5acf0a77-9e40-4117-83fa-4a0791849265
begin
	shift_unit_sel = @bind shift_unit PlutoUI.Slider(sort(unique(spikes.unit)), show_value=true)
	 shift_shift_sel = @bind shift_shift PlutoUI.Slider(sort(collect(keys(shifted))), show_value=true,default=0)
	(;shift_unit_sel, shift_shift_sel)
end

# ╔═╡ 94930aab-8bb0-4da0-b26b-35ddb3efde3b
begin
    obj = Timeshift.ShiftedFields(shifted)
    plot(obj[shift_shift, shift_unit])
end

# ╔═╡ Cell order:
# ╟─0be7ba01-a316-41ce-8df3-a5ae028c74e7
# ╟─37d7f4fd-80a7-47d0-8912-7f002620109f
# ╟─44dde9e4-f9ca-11ec-1348-d968780f671c
# ╟─d33ca749-25af-4b1e-b78c-00ebf54f1327
# ╟─31082fe7-ed61-4d37-a025-77420da3f24a
# ╠═d51ce0f2-03bf-4c88-9302-9fd4bc8621eb
# ╟─efcdc2f1-5e26-4534-953e-defae4bd8603
# ╠═6bdcf863-9946-4ca3-ab02-fa6aebe4b91d
# ╠═03586347-83ee-429d-ab29-505754c66734
# ╟─7c03a7aa-3181-4b3d-9dc6-0eb8e8689023
# ╠═890fe951-19bb-4c9a-a905-d798bb36c57e
# ╠═592d79b4-edf6-4a0c-af73-1d2805d6410e
# ╟─2c794dc4-4920-4274-ab2b-6bb60251112b
# ╟─bef016cd-26d1-4de8-a970-182fe2b92e88
# ╟─410128bf-2332-4aa2-91cb-441e0235cc4a
# ╟─c23ee21a-b3d4-42e5-a4f1-b703008eda1a
# ╠═b88c0ec1-b150-49be-828f-6c32bb770c48
# ╠═38cff24f-bbc1-42fd-98ae-385323c2480e
# ╠═7c6cfeb1-2c78-4480-852b-aa06cc818f76
# ╠═44abcbd4-5f71-4924-b77d-9680cc96044f
# ╟─4d814c3e-97e1-491a-b1d8-c7ca9c628afd
# ╟─f9378d49-2f86-4088-bc6d-3b5b227b7c66
# ╠═5a348e0e-a313-4fe3-a6f7-132fed079c02
# ╠═00736d6f-5a85-4ad7-8661-f2857a36237c
# ╠═f001cdd1-fb4b-4f40-bfcd-be9fea6c8e5b
# ╟─34b5441c-add2-4272-b384-67994daf7745
# ╠═9e5fde4a-f312-419a-b7d5-5ef7a08c305e
# ╠═39737bd9-f38a-408d-a0c0-99b9e2bd0045
# ╟─5acf0a77-9e40-4117-83fa-4a0791849265
# ╠═94930aab-8bb0-4da0-b26b-35ddb3efde3b
