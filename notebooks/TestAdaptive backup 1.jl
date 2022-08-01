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
	  import Plot
	  using Field.metrics
	  
	  adaptive = Field.adaptive
      metrics = Field.metrics
	  WIDTHS = OrderedDict(
		  "x"=>2.5f0, "y"=>2.5f0, "currentPathLength"=>1f0,
          "currentAngle"=>Float32(2pi/80)
	  )
      filts = Filt.get_filters_precache()
	maxrad = nothing
end

# ╔═╡ 5f6e31d3-7101-49aa-a289-39e3967aa3a8
using ColorSchemes

# ╔═╡ cbd7498a-5075-4e67-a829-95f9936146db
using MarkdownLiteral: @mdx

# ╔═╡ 4feb210e-08f1-452a-8d9a-357e12c7210f
using DataFramesMeta

# ╔═╡ e11aae39-41ff-433f-8225-c5230c7c5e2e
using StatsPlots

# ╔═╡ 348e8178-ae24-4217-93a5-54d979b47d92
begin
	using ImageSegmentation, Images, LazySets, Statistics
end

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

# ╔═╡ 04eec95d-b5cf-47a5-a880-3146088cab00
begin
	aspect_ratio = prop_str == "y-x" ? 1 : :none
	ylim = prop_str == "y-x" ? nothing : (0, 100)
	radiusinc = prop_str == "y-x" ? 0.1f0 : 0.05f0
	#thresh = prop_str == "y-x" ? 1.5 : 1
end

# ╔═╡ efcdc2f1-5e26-4534-953e-defae4bd8603
md"""
# Δt = 0 only

## 🌐 Grid
Let's try making a grid object
"""

# ╔═╡ 92b1c56a-1738-43ba-96c9-8c70c6713c39
grid_select

# ╔═╡ 7c03a7aa-3181-4b3d-9dc6-0eb8e8689023
md"""
## 💠 Occupancy
And an occupancy object
"""

# ╔═╡ 2c794dc4-4920-4274-ab2b-6bb60251112b
md"""
## Multiunit adaptive field
grab all spikes as single multiunit field
"""

# ╔═╡ fca06c75-a137-411c-9bd8-74d33ad93183
grid_select

# ╔═╡ c23ee21a-b3d4-42e5-a4f1-b703008eda1a
md"""
## 🔥 Field dict of all cells
"""

# ╔═╡ 38cff24f-bbc1-42fd-98ae-385323c2480e
grid_select

# ╔═╡ 7b150cd8-ada8-46bc-b3ea-1f10c1aea9d8
md"""
### Visualize fields @ Δt=0
"""

# ╔═╡ 7c6cfeb1-2c78-4480-852b-aa06cc818f76
unit_select = @bind unit PlutoUI.Slider(sort(unique(spikes.unit)), default=31, show_value=true)

# ╔═╡ f9378d49-2f86-4088-bc6d-3b5b227b7c66
md"""
### 📜 `to_dataframe`

we're going to want codes that take a set of receptive fields and turn them into a dataframe

- fields
- metrics
"""

# ╔═╡ 52854f4d-1cb2-4a3f-8d48-91aea9c3c45a


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
# ╠═╡ disabled = true
#=╠═╡
I = Timeshift.load_fields();
  ╠═╡ =#

# ╔═╡ b2c8eeb3-9ef7-45ea-926f-26851c7088a4
md"Select a key from prior data?"

# ╔═╡ cbddd54f-83a6-433a-8e0f-67116b476e2e
md"""
New method of picking keys
"""

# ╔═╡ afcb55f9-ed2d-4860-a997-c272bece208f
md"grab that key, or if no key, compute with our settings above"

# ╔═╡ be6048d8-6f30-4d48-a755-5537c3b0b104
md"""
## Visualize timeshifted fields
"""

# ╔═╡ 5731b325-0d82-4141-9122-3b67650ca2a5
md"Select out a unit and a shift"

# ╔═╡ 2a213e93-eaca-416d-9771-788e115e4081
md"Pull out a single shifted field"

# ╔═╡ 15e359aa-6ce2-4479-87cc-4ec3c7a5dfa2
md"plot out a unit at a shift"

# ╔═╡ 47af1633-99bd-4dc2-9d91-9073ec327f27
md"""
## ShiftedField and ShiftedFields objects
"""

# ╔═╡ 8e872986-f947-4c0e-87bb-6e9db34a2508
md"""
Add spatial coherence to this
"""

# ╔═╡ ccb5cb63-4e28-40f5-9866-efceaaf73922


# ╔═╡ 588bff56-6518-4410-b19a-dc745cf067e7
md"""
# Single cell 🦠 metric Development

## Convex Hull 🔴

Want to ensure that watershed segmentation followed by hull of the largest thresholded segments give a good hull
"""

# ╔═╡ f02a79c9-01b8-4550-b321-7b5a6f0d5a28
md"""
### Examine hull creation process
"""

# ╔═╡ 0941a2f5-047f-4a30-823f-fafc53f18b38
segmentation_thresh = @bind qthresh PlutoUI.Slider(0.05:0.05:1; show_value=true, default=0.85)

# ╔═╡ 35da24ff-42c3-455e-9a43-7d74d01f3265
segmentation_thresh

# ╔═╡ d67ba5d1-ca84-4a8b-97ec-54213f60a092
md"""
### Parcellate hull data into metrics
"""

# ╔═╡ 97f2daad-1190-46a2-8c1a-288e0177a29b
function resort_zones(hullzones, instruction::OrderedDict)
	newhullzones = copy(hullzones)
	segsizes = Vector{Int32}(undef, length(instruction))
	for (new, (curr, count)) in enumerate(instruction)
		#if current != new
		#	@info "different"
		#end
		hzinds = hullzones .== Int64(curr)
		#@info curr any(hzinds)
		segsizes[new]         = sum(convert(Array{Int32}, hzinds))
		newhullzones[hzinds] .= new
	end
	return newhullzones, 
		   OrderedDict(zip(1:length(instruction), 
					   values(instruction))),
		   OrderedDict(zip(1:length(instruction),
						   segsizes))
end

# ╔═╡ 87a82eb9-cd22-47f8-acf6-317a794d70ea
segmentation_thresh

# ╔═╡ ce3c4a7c-47d1-4564-bc5c-4c4b20c0820a


# ╔═╡ 3d525c52-6615-4a00-beb6-8c3228f32c6c
md"""
### Hull set operations
Testing about the ability to compute with hulls -- can I see if something is inside a place field
"""

# ╔═╡ 9b9ca9f0-6cce-4f5c-8194-9cd9a2262aa6
md"""
We can either test a set of points for subset relationship
"""

# ╔═╡ 865e9584-6dd9-4d98-bd86-cd918d23379c
md"""
Or single points for set membership
"""

# ╔═╡ c558fd36-5377-4a16-a174-62b8959c4fdc
begin
	struct HullSet
		hulls::Dict
	end
	function Base.:∈(H::HullSet, x)
		any([x ∈ v for v in values(H)])
	end
	function Base.:⊆(H::HullSet, x)
		any([x ⊆ v for v in values(H)])
	end
	Base.iterate(H::HullSet) = Base.iterate(H.hulls)
	Base.iterate(H::HullSet, i::Int64) = Base.iterate(H.hulls, i::Int64)
	Base.keys(H::HullSet) = Base.keys(H.hulls)
	function plothullset!(H::HullSet)
		for i in sort(collect(filter(x-> x isa Int, keys(H))))
			plot!(VPolygon(H.hulls[i]))
			loc = Iterators.flatten(mean(H.hulls[i],dims=1))
			@info loc
			annotate!(loc..., 
					  text(string(i), :white))
		end
	end
end

# ╔═╡ ff355ad4-42da-4493-ae56-3bc9f0d8627c
begin
    widths = OrderedDict(zip(keys(WIDTHS), values(WIDTHS).*width))
    md"""
    widths = $widths
    """
end

# ╔═╡ 6bdcf863-9946-4ca3-ab02-fa6aebe4b91d
# ╠═╡ show_logs = false
# ╠═╡ disabled = true
#=╠═╡
@time G = adaptive.get_grid(beh, props; widths, thresh, maxrad, radiusinc);
  ╠═╡ =#

# ╔═╡ 9e635078-bfdb-41bf-8730-e08a968d5e71
#=╠═╡
md"""
Implied linear width of maxrad[$(nanmaximum(G.radii))] => $(nanmaximum(G.radii) * sqrt(2)) 
"""
  ╠═╡ =#

# ╔═╡ 03586347-83ee-429d-ab29-505754c66734
#=╠═╡
plot(
	plot(G; aspect_ratio, title="radii\nresolution=$(size(G.grid))", ylims=ylim), 
	heatmap([collect(x) for x in G.centers]..., (G.radii .=== NaN32)'; aspect_ratio, title="nan locations")
	, ylims=ylim)
  ╠═╡ =#

# ╔═╡ bf5ec1fc-0443-49df-b90a-164bdd4e8b1b
#=╠═╡
G.centers
  ╠═╡ =#

# ╔═╡ 93b3a5c7-6c6f-4e80-91e7-01f83d292c9a
#=╠═╡
G.radii
  ╠═╡ =#

# ╔═╡ 5550e97c-a33e-4ae4-b888-90f782506bc2
#=╠═╡
unique(G.radii)
  ╠═╡ =#

# ╔═╡ 890fe951-19bb-4c9a-a905-d798bb36c57e
#=╠═╡
O = @time adaptive.get_occupancy(beh, G);
  ╠═╡ =#

# ╔═╡ 592d79b4-edf6-4a0c-af73-1d2805d6410e
#=╠═╡
plot(O, clim=(0,0.01), ylims=ylim)
  ╠═╡ =#

# ╔═╡ d7175827-7528-4cfe-bf3f-d9971f682f49
#=╠═╡
O
  ╠═╡ =#

# ╔═╡ bef016cd-26d1-4de8-a970-182fe2b92e88
#=╠═╡
# Test field abilities
@time multiunit = @time adaptive.get_adaptivefield(spikes, G, O);
# @benchmark adaptive.get_adaptivefield(spikes, G, O);
  ╠═╡ =#

# ╔═╡ 410128bf-2332-4aa2-91cb-441e0235cc4a
#=╠═╡
plot(multiunit)
  ╠═╡ =#

# ╔═╡ b88c0ec1-b150-49be-828f-6c32bb770c48
#=╠═╡
begin
	@time units = adaptive.yartsev(spikes, G, O; widths=width, thresh, 
	                               filters=filts[:all]);
end;
  ╠═╡ =#

# ╔═╡ 44abcbd4-5f71-4924-b77d-9680cc96044f
#=╠═╡
plot(units[(;unit=unit)]; aspect_ratio, ylims=ylim)
  ╠═╡ =#

# ╔═╡ 4d814c3e-97e1-491a-b1d8-c7ca9c628afd
#=╠═╡
μ_firing = begin
    Q = units[(;unit=unit)]
    nansum(reshape(Q.occ.prob, size(Q.occ.count)) .* Q.rate)
end
  ╠═╡ =#

# ╔═╡ eefa56cc-f303-40ee-aa44-dc758eac750b
#=╠═╡
field = units[(;unit=unit)]
  ╠═╡ =#

# ╔═╡ 9405b2bd-c10c-4ba7-aeda-b9f56e2b33ee
# ╠═╡ disabled = true
#=╠═╡
begin
	halfmast = nanquantile(vec(field.rate), qthresh)
	bw = field.rate .> halfmast
	dist = 1 .- distance_transform(feature_transform(bw))
	markers = label_components( (!).(dist .< 0))
end;
  ╠═╡ =#

# ╔═╡ ce81a2d1-7ba8-44fb-b401-760411421a71
#=╠═╡
segments = watershed(dist, markers)
  ╠═╡ =#

# ╔═╡ 8bfb8c40-942d-41b8-a441-5a71f6bbafb7
#=╠═╡
sortperm(collect(values(segments.segment_pixel_count)))


  ╠═╡ =#

# ╔═╡ 794aae46-914a-4da3-a093-d76f1308c55b
#=╠═╡
hullzones = bw .* labels_map(segments);
  ╠═╡ =#

# ╔═╡ 3ce5d298-62eb-4c78-93d8-aa12671fbdce
#=╠═╡
hullzones
  ╠═╡ =#

# ╔═╡ b5e2ef6a-942a-4782-9c36-cd2778de2c66
#=╠═╡
unique(hullzones)
  ╠═╡ =#

# ╔═╡ 99c12e94-8d3e-4700-ab50-146165f654bd
#=╠═╡
plot(
	plot(field, 	 title="field"), 
	heatmap(bw', 	 title="thresholded"), 
	heatmap(dist', 	 title="distance computation"),
	heatmap(markers',title="markers"),
	aspect_ratio=1)
  ╠═╡ =#

# ╔═╡ 0f80805a-76c8-4d8d-91bc-c013575d3a10
#=╠═╡
heatmap(plot(field, title="field"), heatmap(Int8.(hullzones)', title="segments"), aspect_ratio=1)
  ╠═╡ =#


# ╔═╡ d48e9bc9-4d67-4e7b-a0f7-25b975813ccd
#=╠═╡
begin
	
	mets = Dict()

	# Instruction for how to resort by number of pixels
	instruction = OrderedDict(sort([k=>sum(hullzones .== k) for k in 									1:maximum(hullzones)], by=x->x[2], rev=true))
	
	newhullzones, resorted, segsizes = resort_zones(hullzones, instruction)
	
	mets[:hullzone]      = newhullzones
	mets[:hullsegsizes]  = segsizes
	
	ordered_seg = sort(collect(keys(segments.segment_pixel_count)))
	mets[:hullsegsizes] = OrderedDict(k=>segments.segment_pixel_count[k]
									  for k in ordered_seg)
		
	array_of_tuples(X::Vector{<:CartesianIndex})    = [x.I for x in X]
	array_of_arrays(X::Vector{<:CartesianIndex})    = [collect(x.I) for x in X]
	array_of_singleton(X::Vector{<:CartesianIndex}) = [Singleton(collect(x.I)) for x in X]
	
	loopup_coord(c::Tuple, F::Field.ReceptiveField) = F.grid.grid[c...]
	to_grid(X::Vector{<:Union{Tuple,Vector}}, grid::T 
		where T <: Field.adaptive.GridAdaptive) = [grid.grid[Int32.(x)...] for x in X]

	function hull_withimages(X::BitMatrix)::Vector{CartesianIndex}
		convexhull(X)
	end
	function hull_withlazysets(X::BitMatrix)::Vector{Vector{Float32}}
		collect(convex_hull(array_of_arrays(findall(X))))
		#[convert(Vector{Float32}, x) for x in h]
	end
	function centroid(X::BitArray)::Vector{Int32}
		round.(mean(array_of_arrays(findall(X))))
	end
	function centroid(X::BitArray, grid::Field.adaptive.GridAdaptive)::Vector{Float32}
		grid.grid[ centroid(X)... ]
	end

	zones = 1:maximum(newhullzones)
	mets[:hullseg_inds] = Dict{Union{Int, Symbol},Any}()
	mets[:hullseg_grid] = Dict{Union{Int, Symbol},Any}()
	mets[:hullseg_inds_cent] = Dict{Union{Int, Symbol},Any}()
	mets[:hullseg_grid_cent] = Dict{Union{Int, Symbol},Any}()
		
	for zone in zones
		iszone = newhullzones .== zone
		mets[:hullseg_inds][zone] = hull_withlazysets(iszone)
		mets[:hullseg_grid][zone] = to_grid(mets[:hullseg_inds][zone], 
			field.grid)
		mets[:hullseg_inds_cent][zone] = centroid(iszone)
		mets[:hullseg_grid_cent][zone] = centroid(iszone, field.grid)
	end

	iszone = newhullzones .== ordered_seg[1] .||
								newhullzones .== ordered_seg[2]
	mets[:hullseg_inds][:toptwohull]  = hull_withlazysets(iszone)
	mets[:hullseg_grid][:toptwohull] = to_grid(mets[:hullseg_inds][:toptwohull], 												field.grid)
	
	mets[:hullseg_grid_cent][:toptwohull] = centroid(iszone, field.grid)
	mets[:hullseg_inds_cent][:toptwohull] = centroid(iszone)
	mets
	
end;

  ╠═╡ =#

# ╔═╡ 333fe5bf-4168-4931-8856-987d7e76e265
#=╠═╡
instruction
  ╠═╡ =#

# ╔═╡ efd65ca8-457c-4f83-9aaf-162c089404c5
#=╠═╡
plot(heatmap(newhullzones'),heatmap(hullzones'))
  ╠═╡ =#

# ╔═╡ 346a2e86-47be-47ca-9895-fbf4806fc17a
#=╠═╡
ordered_seg
  ╠═╡ =#

# ╔═╡ 70d62e08-77b4-407d-bad8-4850abf5f00a
#=╠═╡
h = hull_withlazysets(hullzones .== 1)
  ╠═╡ =#

# ╔═╡ 8ce9d392-b7bd-483f-b87a-78c6f7657024
#=╠═╡
typeof(h), typeof([h...])
  ╠═╡ =#

# ╔═╡ 0dd5dfe7-321d-466e-9799-ac6c40cc8fb0
#=╠═╡
begin
	plot(VPolygon(h))
	plot!([Singleton(hh) for hh in h], markersize=20)
	plot!(Singleton([3.2f0,8.2f0]), markersize=20)
end
  ╠═╡ =#

# ╔═╡ 30af0459-297b-4c57-95f9-24436d57209c
#=╠═╡
Singleton([3.2f0,8.2f0]) ⊇ VPolygon(h)
  ╠═╡ =#

# ╔═╡ b6b7f2e0-68f7-4324-a78b-197b4143339c
#=╠═╡
Singleton([3.2f0,8.2f0]) ⊆ VPolygon(h)
  ╠═╡ =#

# ╔═╡ 05d6904a-1dfa-4a2f-a385-aaafacc80b0a
#=╠═╡
element(Singleton([3.2f0,8.2f0])) ∈ VPolygon(h)
  ╠═╡ =#

# ╔═╡ d5569288-bc83-408f-8219-5d945cbc6871
segmentation_thresh

# ╔═╡ 0030f529-dd82-4269-aa50-02cc832b9f07
#=╠═╡
begin
	p_with_seghulls = plot(field, aspect_ratio=1)
	plothullset!(HullSet(mets[:hullseg_grid]))
	p_with_seghulls
end
  ╠═╡ =#

# ╔═╡ e06ae752-f36b-465c-9303-74d406b915bf
#=╠═╡
begin
	p_with_seghulls_top = plot(field, aspect_ratio=1)
	plot!(VPolygon(mets[:hullseg_grid][:toptwohull]))
	annotate!(mets[:hullseg_grid_cent][:toptwohull]..., text(string(:toptwohull), :white))
	p_with_seghulls_top
end
  ╠═╡ =#

# ╔═╡ b101021d-e065-4593-b39a-3fee7dbbaf83
#=╠═╡
element(Singleton([50,125])) ∈ HullSet(mets[:hullseg_grid])
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─ff1db172-c3ab-41ea-920c-1dbf831c1336
# ╟─0be7ba01-a316-41ce-8df3-a5ae028c74e7
# ╟─37d7f4fd-80a7-47d0-8912-7f002620109f
# ╠═44dde9e4-f9ca-11ec-1348-d968780f671c
# ╠═5f6e31d3-7101-49aa-a289-39e3967aa3a8
# ╟─a1ad6173-5ead-4559-bddb-9aee6119b9d0
# ╟─2f8ac703-417c-4360-a619-e799d8bb594f
# ╟─31082fe7-ed61-4d37-a025-77420da3f24a
# ╟─d51ce0f2-03bf-4c88-9302-9fd4bc8621eb
# ╟─ff355ad4-42da-4493-ae56-3bc9f0d8627c
# ╟─04eec95d-b5cf-47a5-a880-3146088cab00
# ╟─efcdc2f1-5e26-4534-953e-defae4bd8603
# ╠═6bdcf863-9946-4ca3-ab02-fa6aebe4b91d
# ╟─9e635078-bfdb-41bf-8730-e08a968d5e71
# ╟─92b1c56a-1738-43ba-96c9-8c70c6713c39
# ╠═03586347-83ee-429d-ab29-505754c66734
# ╠═bf5ec1fc-0443-49df-b90a-164bdd4e8b1b
# ╠═93b3a5c7-6c6f-4e80-91e7-01f83d292c9a
# ╠═5550e97c-a33e-4ae4-b888-90f782506bc2
# ╟─7c03a7aa-3181-4b3d-9dc6-0eb8e8689023
# ╠═890fe951-19bb-4c9a-a905-d798bb36c57e
# ╠═592d79b4-edf6-4a0c-af73-1d2805d6410e
# ╠═d7175827-7528-4cfe-bf3f-d9971f682f49
# ╟─2c794dc4-4920-4274-ab2b-6bb60251112b
# ╠═bef016cd-26d1-4de8-a970-182fe2b92e88
# ╟─fca06c75-a137-411c-9bd8-74d33ad93183
# ╟─410128bf-2332-4aa2-91cb-441e0235cc4a
# ╠═c23ee21a-b3d4-42e5-a4f1-b703008eda1a
# ╠═b88c0ec1-b150-49be-828f-6c32bb770c48
# ╟─38cff24f-bbc1-42fd-98ae-385323c2480e
# ╟─7b150cd8-ada8-46bc-b3ea-1f10c1aea9d8
# ╟─7c6cfeb1-2c78-4480-852b-aa06cc818f76
# ╟─44abcbd4-5f71-4924-b77d-9680cc96044f
# ╟─4d814c3e-97e1-491a-b1d8-c7ca9c628afd
# ╠═eefa56cc-f303-40ee-aa44-dc758eac750b
# ╟─f9378d49-2f86-4088-bc6d-3b5b227b7c66
# ╟─52854f4d-1cb2-4a3f-8d48-91aea9c3c45a
# ╟─34b5441c-add2-4272-b384-67994daf7745
# ╟─b2436c80-7290-416f-87c9-137cfb601588
# ╠═cbd7498a-5075-4e67-a829-95f9936146db
# ╠═cb4d5494-24a6-4dfc-980b-23ec48fca7cc
# ╟─b2c8eeb3-9ef7-45ea-926f-26851c7088a4
# ╠═6b8666fc-1e1f-48bf-88bf-b51eb07ad3ce
# ╟─cbddd54f-83a6-433a-8e0f-67116b476e2e
# ╟─c20dd185-4a34-433d-9775-f88475514add
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
# ╠═94930aab-8bb0-4da0-b26b-35ddb3efde3b
# ╟─47af1633-99bd-4dc2-9d91-9073ec327f27
# ╠═62be0c18-ba11-40ad-ba1b-86ec2f638cf3
# ╠═ac9cddcb-097c-42a8-bd59-19fced23bf5a
# ╠═e11aae39-41ff-433f-8225-c5230c7c5e2e
# ╠═e070f95c-1671-4fd7-85b4-20b773de915e
# ╠═8e872986-f947-4c0e-87bb-6e9db34a2508
# ╠═4e394823-6460-4aa7-af24-825a3f379a58
# ╠═d0be8e9b-5f26-46b0-96a0-7ba140b3135f
# ╠═123070a2-0815-49f5-855a-2b87a120e16b
# ╠═ccb5cb63-4e28-40f5-9866-efceaaf73922
# ╟─588bff56-6518-4410-b19a-dc745cf067e7
# ╠═348e8178-ae24-4217-93a5-54d979b47d92
# ╟─f02a79c9-01b8-4550-b321-7b5a6f0d5a28
# ╠═9405b2bd-c10c-4ba7-aeda-b9f56e2b33ee
# ╠═0941a2f5-047f-4a30-823f-fafc53f18b38
# ╠═99c12e94-8d3e-4700-ab50-146165f654bd
# ╠═ce81a2d1-7ba8-44fb-b401-760411421a71
# ╠═794aae46-914a-4da3-a093-d76f1308c55b
# ╠═35da24ff-42c3-455e-9a43-7d74d01f3265
# ╠═0f80805a-76c8-4d8d-91bc-c013575d3a10
# ╟─8bfb8c40-942d-41b8-a441-5a71f6bbafb7
# ╟─d67ba5d1-ca84-4a8b-97ec-54213f60a092
# ╠═3ce5d298-62eb-4c78-93d8-aa12671fbdce
# ╠═333fe5bf-4168-4931-8856-987d7e76e265
# ╠═97f2daad-1190-46a2-8c1a-288e0177a29b
# ╠═d48e9bc9-4d67-4e7b-a0f7-25b975813ccd
# ╠═87a82eb9-cd22-47f8-acf6-317a794d70ea
# ╠═efd65ca8-457c-4f83-9aaf-162c089404c5
# ╠═346a2e86-47be-47ca-9895-fbf4806fc17a
# ╠═b5e2ef6a-942a-4782-9c36-cd2778de2c66
# ╠═70d62e08-77b4-407d-bad8-4850abf5f00a
# ╠═8ce9d392-b7bd-483f-b87a-78c6f7657024
# ╟─ce3c4a7c-47d1-4564-bc5c-4c4b20c0820a
# ╟─3d525c52-6615-4a00-beb6-8c3228f32c6c
# ╟─0dd5dfe7-321d-466e-9799-ac6c40cc8fb0
# ╟─9b9ca9f0-6cce-4f5c-8194-9cd9a2262aa6
# ╠═b6b7f2e0-68f7-4324-a78b-197b4143339c
# ╠═30af0459-297b-4c57-95f9-24436d57209c
# ╟─865e9584-6dd9-4d98-bd86-cd918d23379c
# ╠═05d6904a-1dfa-4a2f-a385-aaafacc80b0a
# ╠═c558fd36-5377-4a16-a174-62b8959c4fdc
# ╠═d5569288-bc83-408f-8219-5d945cbc6871
# ╟─0030f529-dd82-4269-aa50-02cc832b9f07
# ╠═e06ae752-f36b-465c-9303-74d406b915bf
# ╠═b101021d-e065-4593-b39a-3fee7dbbaf83
