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
	  using ColorSchemes
      using Statistics

	  using GoalFetchAnalysis
	  import Utils
	  import Timeshift
	  import Plot
	  using Field.metrics
	  
	  adaptive = Field.adaptive
      metrics  = Field.metrics
	  WIDTHS = OrderedDict(
		  "x"=>2.5f0, "y"=>2.5f0, "currentPathLength"=>1f0,
          "currentAngle"=>Float32(2pi/80)
	  )
      filts = Filt.get_filters_precache()
	  maxrad = nothing

end

# ╔═╡ cbd7498a-5075-4e67-a829-95f9936146db
using MarkdownLiteral: @mdx

# ╔═╡ 4feb210e-08f1-452a-8d9a-357e12c7210f
using DataFramesMeta

# ╔═╡ e11aae39-41ff-433f-8225-c5230c7c5e2e
using StatsPlots

# ╔═╡ ff1db172-c3ab-41ea-920c-1dbf831c1336
md"""
#### 🚀 **Adaptive receptive fields**

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
prop_sel = @bind prop_str PlutoUI.Radio(["y-x","currentAngle-currentPathLength"], default="y-x")

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
	width_select =  @bind width  Slider(0f0:0.2f0:2f0, show_value=true, default=2f0)
	thresh_select = @bind thresh Slider(1f0:1f0:6f0, show_value=true, default=1.5f0)
	(;width_select, thresh_select)
end

# ╔═╡ ff355ad4-42da-4493-ae56-3bc9f0d8627c
begin
    widths = OrderedDict(zip(keys(WIDTHS), values(WIDTHS).*width))
    md"""
    widths = $widths
    """
end

# ╔═╡ 04eec95d-b5cf-47a5-a880-3146088cab00
begin
	aspect_ratio = prop_str == "y-x" ? 1 : :none
	ylim = prop_str == "y-x" ? nothing : (0, 100)
	radiusinc = prop_str == "y-x" ? 0.1f0 : 0.05f0
	#thresh = prop_str == "y-x" ? 1.5 : 1
end

# ╔═╡ 34b5441c-add2-4272-b384-67994daf7745
md"""
# 🔥🕘🕑 Timeshift
Getting everything working with Timeshift.jl
"""

# ╔═╡ b2436c80-7290-416f-87c9-137cfb601588
md"""
First, we run the shifted field calculation
"""

# ╔═╡ cb4d5494-24a6-4dfc-980b-23ec48fca7cc
F = Timeshift.load_fields();

# ╔═╡ b2c8eeb3-9ef7-45ea-926f-26851c7088a4
md"Select a key from prior data?"

# ╔═╡ 6b8666fc-1e1f-48bf-88bf-b51eb07ad3ce
# ╠═╡ disabled = true
#=╠═╡
# OLD METHOD OF PICKING KEYS
begin
	keysets = string.(collect(filter(k->k.grid .== :adaptive .&& k.first.==-2.0 && :Widths ∉ propertynames(k), keys(I))))
	dataset_pick = @bind k PlutoUI.Radio(keysets, default=keysets[2])
end;
  ╠═╡ =#

# ╔═╡ cbddd54f-83a6-433a-8e0f-67116b476e2e
md"""
New method of picking keys
"""

# ╔═╡ c20dd185-4a34-433d-9775-f88475514add
begin
	allsets   = collect(keys(F))
	totalkeys = union(keys.(allsets)...)
	uvals = Dict()
	for key in totalkeys
		push!(uvals, key=>unique([getindex(s, key) for s in allsets if key ∈ propertynames(s)]))
	end
	actuals  = Dict()
	controls = Dict()
	rep(x) = replace(x, "\""=>"", "]"=>"", "["=>"")
	for (i,(K,V)) in enumerate(uvals)
		push!(controls, K => [rep("$v") for v in sort(V)])
		for v in sort(V)
			push!(actuals,  "$v" => v)
		end
	end
	hard_settings = (;grid=:adaptive, first=-2.0, step=0.05, last=2.0)
end

# ╔═╡ 14daafe3-9df8-491d-a7ec-004fa7761c6b
controls

# ╔═╡ c461e634-d011-49bb-9ba5-b01df547f36f
begin
	datacut_sel = @bind datacut_key PlutoUI.Radio(controls[:datacut], default="all")
	thresh_sel  = @bind thresh_key PlutoUI.Radio(controls[:thresh], default="1.5")
	widths_sel  = @bind widths_key PlutoUI.Radio(controls[:widths], default="5.0")
	console = (;datacut_sel, thresh_sel, widths_sel)
end

# ╔═╡ d423498d-d024-4465-8dae-9eb899a75457
key = Utils.namedtup.bestpartialmatch(keys(F), 
		(;hard_settings..., datacut=actuals[datacut_key], thresh=actuals[thresh_key], widths=actuals[widths_key]))

# ╔═╡ afcb55f9-ed2d-4860-a997-c272bece208f
md"grab that key, or if no key, compute with our settings above"

# ╔═╡ 955c7b75-00d4-4116-9242-92a7df8a0f87
shifted = F[key];

# ╔═╡ 2331218a-bc76-4adf-82e3-8e5b52aef0ca
plot_obj = Timeshift.DictOfShiftOfUnit{Float64}(shifted);

# ╔═╡ 5f15dc20-cf30-4088-a173-9c084ac2809a
# ╠═╡ show_logs = false
@time SFs = Timeshift.ShiftedFields(plot_obj);

# ╔═╡ 1994b7d7-c152-4aa7-a088-3f272246d9ae
nonzero_units = sort(@subset(SFs.metrics,
                             :metric .== Symbol("maxcount"),
                             getindex.(:value,1) .!= 0,
                             :shift.==  0).unit)

# ╔═╡ be6048d8-6f30-4d48-a755-5537c3b0b104
md"""
## Visualize timeshifted fields
"""

# ╔═╡ 5731b325-0d82-4141-9122-3b67650ca2a5
md"Select out a unit and a shift"

# ╔═╡ 5acf0a77-9e40-4117-83fa-4a0791849265
begin
	unit_sel = @bind shift_unit PlutoUI.Slider(nonzero_units, show_value=true)
	 shift_sel = @bind shift_shift PlutoUI.Slider(sort(collect(keys(shifted))), show_value=true,default=0)
	unit_shift_sel = (;unit_sel, shift_sel)
end

# ╔═╡ 8c91224f-e205-464c-a088-57393bc99d13
sort(collect(keys(shifted)))

# ╔═╡ 2a213e93-eaca-416d-9771-788e115e4081
md"Pull out a single shifted field"

# ╔═╡ 3afa8aae-bf3e-4364-8ad9-76fd50ca5ac9
single_shifted_field = get(plot_obj, shift_shift, shift_unit);

# ╔═╡ 15e359aa-6ce2-4479-87cc-4ec3c7a5dfa2
md"plot out a unit at a shift"

# ╔═╡ 89e8dc9a-b7d9-4d1f-b604-8376b28bfca5
(;console,unit_shift_sel)

# ╔═╡ dc714977-349e-417f-a5fd-1bd3a8c1e38b
md"""
key = $key
"""

# ╔═╡ 534bca46-d351-4385-97b9-f0358c342fb6
item = get(plot_obj,shift_shift, shift_unit);

# ╔═╡ 94930aab-8bb0-4da0-b26b-35ddb3efde3b
begin
    plot(item; 
		aspect_ratio, ylims=ylim,
		title=string(item.metrics)
	)
end

# ╔═╡ 47af1633-99bd-4dc2-9d91-9073ec327f27
md"""
## ShiftedField and ShiftedFields objects
"""

# ╔═╡ 62be0c18-ba11-40ad-ba1b-86ec2f638cf3
SF = SFs[shift_unit];

# ╔═╡ 16869357-2e3d-4826-b408-a99d5cdf8492
SF.metrics

# ╔═╡ ac9cddcb-097c-42a8-bd59-19fced23bf5a
begin
    shmet = sort(unstack(SF.metrics, :shift, :metric, :value, allowduplicates=true),:shift)
    Table.vec_arrayofarrays!(shmet)
end

# ╔═╡ e070f95c-1671-4fd7-85b4-20b773de915e
(;console,unit_shift_sel)

# ╔═╡ 8e872986-f947-4c0e-87bb-6e9db34a2508
md"""
Add spatial coherence to this
"""

# ╔═╡ 4e394823-6460-4aa7-af24-825a3f379a58
if :information ∈ propertynames(shmet)
	@df shmet histogram(:information)
else
	
	pbps_hist = @df shmet histogram(:bitsperspike, xlabel="bits/spike", label="")
	pbps_shift_scat = @df  shmet scatter(:shift, :bitsperspike, label="", c=get(ColorSchemes.acton, Utils.norm_extrema(:shift,[0,1])))
	shiftopt = shmet.shift[argmax(shmet.bitsperspike)]
	pbps_max = vline!([shiftopt], c=:gray, label="best")
	vline!([0], c=:black, linestyle=:dash, label="")
	annotate!(shiftopt, mean(ylims()), text(shiftopt))

	pc_hist = @df shmet histogram(:coherence, label="", xlabel="spatial coherence")
	pc_shift_scat = @df  shmet scatter(:shift, :coherence, label="", c=get(ColorSchemes.acton, Utils.norm_extrema(:shift,[0,1])))
	shiftopt = shmet.shift[argmax(shmet.coherence)]
	vline!([0], c=:black, linestyle=:dash, label="")
	vline!([shiftopt], c=:gray, label="best")
	annotate!(shiftopt, mean(ylims()), text(shiftopt))
	
	
	plot(pbps_hist, pbps_shift_scat, pc_hist, pc_shift_scat)
end


# ╔═╡ d0be8e9b-5f26-46b0-96a0-7ba140b3135f
@df shmet scatter(:coherence, :bitsperspike, c=get(ColorSchemes.acton, Utils.norm_extrema(:shift,[0,1])), xlabel="Spatial Coherence", ylabel="Bits per spike", label="Unit = $shift_unit", legend=:topleft)

# ╔═╡ 123070a2-0815-49f5-855a-2b87a120e16b
@df shmet scatter3d(:shift, :coherence, :bitsperspike, c=get(ColorSchemes.acton, Utils.norm_extrema(:shift,[0,1])), xlabel="shift",ylabel="coherence", zlabel="bitsperspike")

# ╔═╡ 348e8178-ae24-4217-93a5-54d979b47d92


# ╔═╡ Cell order:
# ╟─ff1db172-c3ab-41ea-920c-1dbf831c1336
# ╟─0be7ba01-a316-41ce-8df3-a5ae028c74e7
# ╟─37d7f4fd-80a7-47d0-8912-7f002620109f
# ╠═44dde9e4-f9ca-11ec-1348-d968780f671c
# ╟─a1ad6173-5ead-4559-bddb-9aee6119b9d0
# ╟─2f8ac703-417c-4360-a619-e799d8bb594f
# ╟─31082fe7-ed61-4d37-a025-77420da3f24a
# ╟─d51ce0f2-03bf-4c88-9302-9fd4bc8621eb
# ╟─ff355ad4-42da-4493-ae56-3bc9f0d8627c
# ╟─04eec95d-b5cf-47a5-a880-3146088cab00
# ╟─34b5441c-add2-4272-b384-67994daf7745
# ╟─b2436c80-7290-416f-87c9-137cfb601588
# ╠═cbd7498a-5075-4e67-a829-95f9936146db
# ╠═cb4d5494-24a6-4dfc-980b-23ec48fca7cc
# ╟─b2c8eeb3-9ef7-45ea-926f-26851c7088a4
# ╠═6b8666fc-1e1f-48bf-88bf-b51eb07ad3ce
# ╟─cbddd54f-83a6-433a-8e0f-67116b476e2e
# ╠═c20dd185-4a34-433d-9775-f88475514add
# ╠═14daafe3-9df8-491d-a7ec-004fa7761c6b
# ╠═c461e634-d011-49bb-9ba5-b01df547f36f
# ╠═d423498d-d024-4465-8dae-9eb899a75457
# ╟─afcb55f9-ed2d-4860-a997-c272bece208f
# ╠═955c7b75-00d4-4116-9242-92a7df8a0f87
# ╠═2331218a-bc76-4adf-82e3-8e5b52aef0ca
# ╠═5f15dc20-cf30-4088-a173-9c084ac2809a
# ╠═4feb210e-08f1-452a-8d9a-357e12c7210f
# ╠═1994b7d7-c152-4aa7-a088-3f272246d9ae
# ╟─be6048d8-6f30-4d48-a755-5537c3b0b104
# ╟─5731b325-0d82-4141-9122-3b67650ca2a5
# ╠═5acf0a77-9e40-4117-83fa-4a0791849265
# ╠═8c91224f-e205-464c-a088-57393bc99d13
# ╟─2a213e93-eaca-416d-9771-788e115e4081
# ╠═3afa8aae-bf3e-4364-8ad9-76fd50ca5ac9
# ╟─15e359aa-6ce2-4479-87cc-4ec3c7a5dfa2
# ╠═89e8dc9a-b7d9-4d1f-b604-8376b28bfca5
# ╟─dc714977-349e-417f-a5fd-1bd3a8c1e38b
# ╠═534bca46-d351-4385-97b9-f0358c342fb6
# ╠═94930aab-8bb0-4da0-b26b-35ddb3efde3b
# ╟─47af1633-99bd-4dc2-9d91-9073ec327f27
# ╠═62be0c18-ba11-40ad-ba1b-86ec2f638cf3
# ╠═16869357-2e3d-4826-b408-a99d5cdf8492
# ╠═ac9cddcb-097c-42a8-bd59-19fced23bf5a
# ╠═e11aae39-41ff-433f-8225-c5230c7c5e2e
# ╠═e070f95c-1671-4fd7-85b4-20b773de915e
# ╠═8e872986-f947-4c0e-87bb-6e9db34a2508
# ╠═4e394823-6460-4aa7-af24-825a3f379a58
# ╠═d0be8e9b-5f26-46b0-96a0-7ba140b3135f
# ╠═123070a2-0815-49f5-855a-2b87a120e16b
# ╠═348e8178-ae24-4217-93a5-54d979b47d92
