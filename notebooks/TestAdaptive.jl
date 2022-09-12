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
# ╠═╡ show_logs = false
begin
	  using DrWatson
	  quickactivate(expanduser("~/Projects/goal-code"))
	  using Plots
	  using Revise
	  using DataFrames
	  using NaNStatistics
	  import ProgressLogging
	  using PlutoUI
	  using DataStructures: OrderedDict

	  using GoalFetchAnalysis
	  import Utils
	  import Timeshift
	  #import Plot
	  using Field.metrics
    using ColorSchemes
    using DataFramesMeta
	using Timeshift.types
	using Timeshift.shiftmetrics
    using StatsBase
    using DimensionalData
    using StatsPlots
    using GLM
	using ImageSegmentation, Images, LazySets, Statistics
    import Plot
	  
	  adaptive = Field.adaptive
      metrics = Field.metrics
	  WIDTHS = OrderedDict(
		  "x"=>2.5f0, "y"=>2.5f0, "currentPathLength"=>2f0,
          "currentAngle"=>Float32(2pi/80)
	  )
      filts = Filt.get_filters_precache()
	maxrad = nothing
end

# ╔═╡ ff1db172-c3ab-41ea-920c-1dbf831c1336
md"""
!!! notebook
	🚀 **Adaptive receptive fields**

Purpose: Test out my base adaptive field codes. Make sure the various steps (grid, occupancy, and fields) are working. Also make sure downstream shifted objects are working.

##### *TODO*
* metrics(shifted fields)
* shifted field plot recipes
* radii vector support
  * plot recipe
  * selection via 10 increase in base sample width
* to_dataframe
  * fields
  * metrics
  * shifted fields
"""

# ╔═╡ 0be7ba01-a316-41ce-8df3-a5ae028c74e7
PlutoUI.TableOfContents(title="🚀 Adaptive RFs" )

# ╔═╡ 37d7f4fd-80a7-47d0-8912-7f002620109f
md"""
# Preamble 
Import packages
"""

# ╔═╡ a1ad6173-5ead-4559-bddb-9aee6119b9d0
prop_sel = @bind prop_str PlutoUI.Radio(["y-x","currentAngle-currentPathLength", "currentAngle","currentPathLength"], default="y-x")

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
	width_select =  @bind width  Slider(0f0:0.2f0:2f0, show_value=true, default=04f0)
	thresh_select = @bind thresh Slider(1f0:0.5f0:6f0, show_value=true, default=4f0)
	(;width_select, thresh_select)
end

# ╔═╡ ff355ad4-42da-4493-ae56-3bc9f0d8627c
begin
    widths = OrderedDict(zip(keys(WIDTHS), values(WIDTHS).*width))
    md"""
    widths = $widths
    """
end

# ╔═╡ cbfbe9c9-f7c5-4219-8d81-b335fe5f5ed6
radiusinc, ylim, aspect_ratio = if prop_sel == "y-x"
	0.1f0, nothing, 1
elseif prop_sel == "currentAngle-currentPathLength"
    0.05f0, (0, 100), 1/18
else
	0.05f0, nothing, 1
end

# ╔═╡ efcdc2f1-5e26-4534-953e-defae4bd8603
md"""
# Δt = 0 only

## 🌐 Grid
Let's try making a grid object
"""

# ╔═╡ 301d385b-7879-445d-a903-772c10862750
widths

# ╔═╡ 6bdcf863-9946-4ca3-ab02-fa6aebe4b91d
# ╠═╡ show_logs = false
begin
    boundary = prop_str == "y-x" ? nothing : Dict("currentPathLength"=>(0,100))
    @time G = adaptive.get_grid(beh, props; widths, thresh, maxrad=prop_str == "y-x" ? 10f0 : nothing, radiusinc, boundary);
end

# ╔═╡ 9e635078-bfdb-41bf-8730-e08a968d5e71
md"""
Implied linear width of `maxrad[$(nanmaximum(G.radii))]` => $(nanmaximum(G.radii) * sqrt(2)) 
"""

# ╔═╡ 92b1c56a-1738-43ba-96c9-8c70c6713c39
grid_select

# ╔═╡ 03586347-83ee-429d-ab29-505754c66734
plot(G; title="radii\nresolution=$(size(G.grid))")

# ╔═╡ 5dda1153-5359-49d5-9baf-b3bebc4627a0
if ndims(G.radii) == 2
	heatmap([collect(x) for x in G.centers]..., (G.radii .=== NaN32)'; title="nan locations")
end

# ╔═╡ bf5ec1fc-0443-49df-b90a-164bdd4e8b1b
G.centers

# ╔═╡ 93b3a5c7-6c6f-4e80-91e7-01f83d292c9a
G.radii

# ╔═╡ 5550e97c-a33e-4ae4-b888-90f782506bc2
unique(G.radii)

# ╔═╡ 7c03a7aa-3181-4b3d-9dc6-0eb8e8689023
md"""
## 💠 Occupancy
And an occupancy object
"""

# ╔═╡ 890fe951-19bb-4c9a-a905-d798bb36c57e
O = @time adaptive.get_occupancy(beh, G);

# ╔═╡ 592d79b4-edf6-4a0c-af73-1d2805d6410e
plot(O, clim=(0,0.01), ylims=ylim)

# ╔═╡ 5dfeb288-08a2-4393-8eb9-978fadcf57a8


# ╔═╡ d7175827-7528-4cfe-bf3f-d9971f682f49
# ╠═╡ disabled = true
#=╠═╡
O
  ╠═╡ =#

# ╔═╡ 2c794dc4-4920-4274-ab2b-6bb60251112b
md"""
## Multiunit adaptive field
grab all spikes as single multiunit field
"""

# ╔═╡ bef016cd-26d1-4de8-a970-182fe2b92e88
# Test field abilities
@time multiunit = @time adaptive.get_adaptivefield(spikes, G, O);
# @benchmark adaptive.get_adaptivefield(spikes, G, O);

# ╔═╡ fca06c75-a137-411c-9bd8-74d33ad93183
grid_select

# ╔═╡ 410128bf-2332-4aa2-91cb-441e0235cc4a
plot(multiunit)

# ╔═╡ c23ee21a-b3d4-42e5-a4f1-b703008eda1a
md"""
## 🔥 Field dict of all cells
"""

# ╔═╡ b88c0ec1-b150-49be-828f-6c32bb770c48
begin
	@time units = adaptive.yartsev(spikes, G, O; widths=width, thresh, 
	                               filters=filts[:all]);
end;

# ╔═╡ 38cff24f-bbc1-42fd-98ae-385323c2480e
grid_select

# ╔═╡ 7b150cd8-ada8-46bc-b3ea-1f10c1aea9d8
md"""
### Visualize fields @ Δt=0
"""

# ╔═╡ 589f89f4-bc38-49c3-8973-e9e52a26647b
revise(Plot.receptivefield)

# ╔═╡ bb23f597-480e-43d1-a527-71f05a426d1d
begin
	bps = OrderedDict(k=>v[:bitsperspike] for (k,v) in units)
	sort!(bps, by=k->bps[k], rev=true)
	units_ordered = [x[1] for x in keys(bps)]
end

# ╔═╡ c34cfe1f-3663-4901-a258-0015cfa74a38
# ╠═╡ disabled = true
#=╠═╡
Memoization.empty_all_caches!()
  ╠═╡ =#

# ╔═╡ a0d00fd1-c587-4bee-aa14-d0366a7f65ae
# ╠═╡ disabled = true
#=╠═╡
Plot.setfolder("goal", "pathlength")
  ╠═╡ =#

# ╔═╡ 7c6cfeb1-2c78-4480-852b-aa06cc818f76
unit_select = @bind unit PlutoUI.Slider(units_ordered, default=31, show_value=true)

# ╔═╡ cc1c449e-2b11-40c0-a24c-5edc1ba69457
U = units[(;unit=unit)];

# ╔═╡ 44abcbd4-5f71-4924-b77d-9680cc96044f
plot(U)

# ╔═╡ 3f6cec15-f43d-43d2-866b-45a04ea5715f
# ╠═╡ disabled = true
#=╠═╡
Plot.save((;desc="pathlength",unit))
  ╠═╡ =#

# ╔═╡ f2fc5e48-3adc-41d0-b3e9-36bca362415b
ylims === nothing

# ╔═╡ 4d814c3e-97e1-491a-b1d8-c7ca9c628afd
μ_firing = begin
    Q = units[(;unit=unit)]
    nansum(reshape(Q.occ.prob, size(Q.occ.count)) .* Q.rate)
end

# ╔═╡ eefa56cc-f303-40ee-aa44-dc758eac750b
field = units[(;unit=unit)]

# ╔═╡ Cell order:
# ╟─ff1db172-c3ab-41ea-920c-1dbf831c1336
# ╟─0be7ba01-a316-41ce-8df3-a5ae028c74e7
# ╟─37d7f4fd-80a7-47d0-8912-7f002620109f
# ╠═44dde9e4-f9ca-11ec-1348-d968780f671c
# ╠═a1ad6173-5ead-4559-bddb-9aee6119b9d0
# ╟─2f8ac703-417c-4360-a619-e799d8bb594f
# ╟─31082fe7-ed61-4d37-a025-77420da3f24a
# ╠═d51ce0f2-03bf-4c88-9302-9fd4bc8621eb
# ╟─ff355ad4-42da-4493-ae56-3bc9f0d8627c
# ╠═cbfbe9c9-f7c5-4219-8d81-b335fe5f5ed6
# ╟─efcdc2f1-5e26-4534-953e-defae4bd8603
# ╟─301d385b-7879-445d-a903-772c10862750
# ╠═6bdcf863-9946-4ca3-ab02-fa6aebe4b91d
# ╟─9e635078-bfdb-41bf-8730-e08a968d5e71
# ╟─92b1c56a-1738-43ba-96c9-8c70c6713c39
# ╠═03586347-83ee-429d-ab29-505754c66734
# ╠═5dda1153-5359-49d5-9baf-b3bebc4627a0
# ╠═bf5ec1fc-0443-49df-b90a-164bdd4e8b1b
# ╠═93b3a5c7-6c6f-4e80-91e7-01f83d292c9a
# ╠═5550e97c-a33e-4ae4-b888-90f782506bc2
# ╟─7c03a7aa-3181-4b3d-9dc6-0eb8e8689023
# ╠═890fe951-19bb-4c9a-a905-d798bb36c57e
# ╠═592d79b4-edf6-4a0c-af73-1d2805d6410e
# ╠═5dfeb288-08a2-4393-8eb9-978fadcf57a8
# ╠═d7175827-7528-4cfe-bf3f-d9971f682f49
# ╟─2c794dc4-4920-4274-ab2b-6bb60251112b
# ╠═bef016cd-26d1-4de8-a970-182fe2b92e88
# ╟─fca06c75-a137-411c-9bd8-74d33ad93183
# ╠═410128bf-2332-4aa2-91cb-441e0235cc4a
# ╟─c23ee21a-b3d4-42e5-a4f1-b703008eda1a
# ╠═b88c0ec1-b150-49be-828f-6c32bb770c48
# ╟─38cff24f-bbc1-42fd-98ae-385323c2480e
# ╟─7b150cd8-ada8-46bc-b3ea-1f10c1aea9d8
# ╠═cc1c449e-2b11-40c0-a24c-5edc1ba69457
# ╠═589f89f4-bc38-49c3-8973-e9e52a26647b
# ╠═bb23f597-480e-43d1-a527-71f05a426d1d
# ╠═c34cfe1f-3663-4901-a258-0015cfa74a38
# ╠═a0d00fd1-c587-4bee-aa14-d0366a7f65ae
# ╟─7c6cfeb1-2c78-4480-852b-aa06cc818f76
# ╠═44abcbd4-5f71-4924-b77d-9680cc96044f
# ╠═3f6cec15-f43d-43d2-866b-45a04ea5715f
# ╠═f2fc5e48-3adc-41d0-b3e9-36bca362415b
# ╟─4d814c3e-97e1-491a-b1d8-c7ca9c628afd
# ╠═eefa56cc-f303-40ee-aa44-dc758eac750b
