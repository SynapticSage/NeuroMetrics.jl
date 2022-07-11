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
	  using NaNStatistics
	  import ProgressLogging
	  import Utils
	  import Timeshift
	  using PlutoUI
	  using DataStructures: OrderedDict
	  
	  adaptive = Field.adaptive
      metrics = Field.metrics
	  WIDTHS = OrderedDict(
		  "x"=>2.5f0, "y"=>2.5f0, "currentPathLength"=>3f0, "currentAngle"=>Float32(2pi/20)
	  )
end

# ╔═╡ 0be7ba01-a316-41ce-8df3-a5ae028c74e7
PlutoUI.TableOfContents(title="Test adaptive")

# ╔═╡ 37d7f4fd-80a7-47d0-8912-7f002620109f
md"""
# Preamble 
Import packages
"""

# ╔═╡ 165e0caa-b491-43d2-b468-3456c660dd02
valtype(WIDTHS)

# ╔═╡ a1ad6173-5ead-4559-bddb-9aee6119b9d0
prop_sel = @bind prop_str PlutoUI.Radio(["x-y","currentAngle-currentPathLength"], 										default="x-y")

# ╔═╡ 31082fe7-ed61-4d37-a025-77420da3f24a
beh, spikes = begin
	props = Vector{String}(split(prop_str, "-"))
	@info props
	@time spikes, beh, ripples, cells = Load.load("RY16", 36);
	@time beh, spikes  = Load.register(beh, spikes; on="time", transfer=props)
    spikes = dropmissing(spikes, props)
    beh, spikes
end;

# ╔═╡ 2f8ac703-417c-4360-a619-e799d8bb594f
md"""
Loadup (spikes,beh) dataframes and transfer $(join(props,"-")) to spikes structure
"""

# ╔═╡ d51ce0f2-03bf-4c88-9302-9fd4bc8621eb
grid_select = begin
	width_select =  @bind width  Slider(0f0:0.2f0:2f0, show_value=true, default=1f0)
	thresh_select = @bind thresh Slider(1f0:1f0:6f0, show_value=true, default=1.5f0)
	(;width_select, thresh_select)
end

# ╔═╡ ff355ad4-42da-4493-ae56-3bc9f0d8627c
widths = OrderedDict(zip(keys(WIDTHS), values(WIDTHS).*width))

# ╔═╡ efcdc2f1-5e26-4534-953e-defae4bd8603
md"""
# Grid
Let's try making a grid object
"""

# ╔═╡ 6104813f-e8cf-42fe-8f28-16cacf11cce3
maxrad = nothing

# ╔═╡ 6bdcf863-9946-4ca3-ab02-fa6aebe4b91d
G = adaptive.get_grid(beh, props; widths, thresh, maxrad);

# ╔═╡ 9e635078-bfdb-41bf-8730-e08a968d5e71
md"""
Implied linear width of maxrad[$(nanmaximum(G.radii))] => $(nanmaximum(G.radii) * sqrt(2)) 
"""

# ╔═╡ 92b1c56a-1738-43ba-96c9-8c70c6713c39
grid_select

# ╔═╡ 03586347-83ee-429d-ab29-505754c66734
plot(plot(G; aspect_ratio=1, title="radii\nresolution=$(size(G.grid))"), heatmap(isnan.(G.radii), aspect_ratio=1, title="nan locations"))

# ╔═╡ 7c03a7aa-3181-4b3d-9dc6-0eb8e8689023
md"""
# Occupancy
And an occupancy object
"""

# ╔═╡ 890fe951-19bb-4c9a-a905-d798bb36c57e
O = @time adaptive.get_occupancy(beh, G);

# ╔═╡ 592d79b4-edf6-4a0c-af73-1d2805d6410e
plot(O, clim=(0,0.001))

# ╔═╡ 2c794dc4-4920-4274-ab2b-6bb60251112b
md"""
# Single adaptive field
grab all spikes as single multiunit field
"""

# ╔═╡ bef016cd-26d1-4de8-a970-182fe2b92e88
# ╠═╡ disabled = true
#=╠═╡
# Test field abilities
@time multiunit = @time adaptive.get_adaptivefield(spikes, G, O);
  ╠═╡ =#

# ╔═╡ fca06c75-a137-411c-9bd8-74d33ad93183
grid_select

# ╔═╡ 410128bf-2332-4aa2-91cb-441e0235cc4a
#=╠═╡
plot(multiunit)
  ╠═╡ =#

# ╔═╡ c23ee21a-b3d4-42e5-a4f1-b703008eda1a
md"""
# Field dict of all cells
"""

# ╔═╡ b88c0ec1-b150-49be-828f-6c32bb770c48
@time units = adaptive.yartsev(spikes, G, O; splitby=[:unit], widths=width, thresh);

# ╔═╡ 38cff24f-bbc1-42fd-98ae-385323c2480e
grid_select

# ╔═╡ 7b150cd8-ada8-46bc-b3ea-1f10c1aea9d8
md"""
## Visualize fields @ Δt=0
"""

# ╔═╡ 7c6cfeb1-2c78-4480-852b-aa06cc818f76
unit_select = @bind unit PlutoUI.Slider(sort(unique(spikes.unit)), show_value=true)

# ╔═╡ 44abcbd4-5f71-4924-b77d-9680cc96044f
plot(units[(;unit=unit)], aspect_ratio=5)

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

# ╔═╡ b2436c80-7290-416f-87c9-137cfb601588
md"""
First, we run the shifted field calculation
"""

# ╔═╡ 39737bd9-f38a-408d-a0c0-99b9e2bd0045
# ╠═╡ disabled = true
#=╠═╡
shifted = Timeshift.shifted_fields(beh, spikes, -1:1:1, props, widths=width);
  ╠═╡ =#

# ╔═╡ be6048d8-6f30-4d48-a755-5537c3b0b104
md"""
## Visualize timeshifted fields
"""

# ╔═╡ 5acf0a77-9e40-4117-83fa-4a0791849265
#=╠═╡
begin
	unit_sel = @bind shift_unit PlutoUI.Slider(sort(unique(spikes.unit)), show_value=true)
	 shift_sel = @bind shift_shift PlutoUI.Slider(sort(collect(keys(shifted))), show_value=true,default=0)
	(;unit_sel, shift_sel)
end
  ╠═╡ =#

# ╔═╡ 94930aab-8bb0-4da0-b26b-35ddb3efde3b
#=╠═╡
begin
    plot_obj = Timeshift.DictOfShiftOfUnit{Float64}(shifted)
    plot(get(plot_obj,shift_shift, shift_unit), aspect_ratio=1, 
		title=string(get(plot_obj, shift_shift, shift_unit)))
end
  ╠═╡ =#

# ╔═╡ 47af1633-99bd-4dc2-9d91-9073ec327f27
md"""
## ShiftedField and ShiftedFields objects
"""

# ╔═╡ 5f15dc20-cf30-4088-a173-9c084ac2809a
#=╠═╡
S = Timeshift.ShiftedField(get(plot_obj, :, shift_unit))
  ╠═╡ =#

# ╔═╡ dbd958b4-e9c9-4ac4-bb1e-64c5d83f231a


# ╔═╡ Cell order:
# ╠═0be7ba01-a316-41ce-8df3-a5ae028c74e7
# ╟─37d7f4fd-80a7-47d0-8912-7f002620109f
# ╠═44dde9e4-f9ca-11ec-1348-d968780f671c
# ╠═165e0caa-b491-43d2-b468-3456c660dd02
# ╠═a1ad6173-5ead-4559-bddb-9aee6119b9d0
# ╟─2f8ac703-417c-4360-a619-e799d8bb594f
# ╟─31082fe7-ed61-4d37-a025-77420da3f24a
# ╠═d51ce0f2-03bf-4c88-9302-9fd4bc8621eb
# ╠═ff355ad4-42da-4493-ae56-3bc9f0d8627c
# ╟─efcdc2f1-5e26-4534-953e-defae4bd8603
# ╟─6104813f-e8cf-42fe-8f28-16cacf11cce3
# ╠═6bdcf863-9946-4ca3-ab02-fa6aebe4b91d
# ╟─9e635078-bfdb-41bf-8730-e08a968d5e71
# ╟─92b1c56a-1738-43ba-96c9-8c70c6713c39
# ╟─03586347-83ee-429d-ab29-505754c66734
# ╟─7c03a7aa-3181-4b3d-9dc6-0eb8e8689023
# ╟─890fe951-19bb-4c9a-a905-d798bb36c57e
# ╠═592d79b4-edf6-4a0c-af73-1d2805d6410e
# ╟─2c794dc4-4920-4274-ab2b-6bb60251112b
# ╠═bef016cd-26d1-4de8-a970-182fe2b92e88
# ╟─fca06c75-a137-411c-9bd8-74d33ad93183
# ╟─410128bf-2332-4aa2-91cb-441e0235cc4a
# ╟─c23ee21a-b3d4-42e5-a4f1-b703008eda1a
# ╟─b88c0ec1-b150-49be-828f-6c32bb770c48
# ╟─38cff24f-bbc1-42fd-98ae-385323c2480e
# ╟─7b150cd8-ada8-46bc-b3ea-1f10c1aea9d8
# ╟─7c6cfeb1-2c78-4480-852b-aa06cc818f76
# ╠═44abcbd4-5f71-4924-b77d-9680cc96044f
# ╟─4d814c3e-97e1-491a-b1d8-c7ca9c628afd
# ╟─f9378d49-2f86-4088-bc6d-3b5b227b7c66
# ╟─34b5441c-add2-4272-b384-67994daf7745
# ╟─b2436c80-7290-416f-87c9-137cfb601588
# ╟─39737bd9-f38a-408d-a0c0-99b9e2bd0045
# ╟─be6048d8-6f30-4d48-a755-5537c3b0b104
# ╟─5acf0a77-9e40-4117-83fa-4a0791849265
# ╟─94930aab-8bb0-4da0-b26b-35ddb3efde3b
# ╟─47af1633-99bd-4dc2-9d91-9073ec327f27
# ╠═5f15dc20-cf30-4088-a173-9c084ac2809a
# ╠═dbd958b4-e9c9-4ac4-bb1e-64c5d83f231a
