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

# ╔═╡ 5dfeb288-08a2-4393-8eb9-978fadcf57a8


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

# ╔═╡ 589f89f4-bc38-49c3-8973-e9e52a26647b
revise(Plot.receptivefield)

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

# ╔═╡ f2fc5e48-3adc-41d0-b3e9-36bca362415b
ylims === nothing

# ╔═╡ f9378d49-2f86-4088-bc6d-3b5b227b7c66
md"""
### 📜 `to_dataframe`

we're going to want codes that take a set of receptive fields and turn them into a dataframe

- fields
- metrics
"""

# ╔═╡ d1fe3b02-34f5-46ae-8da8-4bac71c86d84
shifts = -2.5:0.2:2.5

# ╔═╡ 014d81aa-5441-41c1-8b71-a347dd03ac9a


# ╔═╡ 316dab79-015f-46e4-88f5-7051529484e5
@bind timeshift Slider(shifts, default=0)

# ╔═╡ 69d37df3-69d3-4c91-8adf-ff8ecd4c9df1
md"## Isolated spiking"

# ╔═╡ 664fba69-e498-462a-8e84-522a22f2498d
if !hasproperty(spikes, :isolated)
	@time lfp = Load.load_lfp("RY16",36,tet=5)
	@time Munge.lfp.annotate_cycles(lfp)
	lfp.time .-= Load.min_time_records[1]
	@time Munge.spiking.isolated(spikes,lfp)
end

# ╔═╡ bf12b73f-c146-4ab2-be79-5b8097add3f3
begin
	datacut_sl = @bind datacut_str PlutoUI.Radio(["all","cue","memory","nontask"], default="memory")

	norm_sl = @bind norm PlutoUI.Radio(["01","percent"],default="percent")
	
	iso_sl=(;datacut_sl, norm_sl)
end

# ╔═╡ 00936498-8061-43da-bdb3-0c7b000770bb
datacut = Symbol(datacut_str)

# ╔═╡ ffe4b420-41ae-4c58-8bd6-d5d684f48b78
if norm == "percent"
		normf = Utils.norm_percent
	else
		normf = Utils.norm_extrema
	end

# ╔═╡ fff27c22-7c2c-45d0-a670-6920b9bd6a31
all((!).(spikes.isolated))

# ╔═╡ 85cf402a-eb7b-4d3c-972b-810ac1708632
isotab = combine(groupby(spikes, :unit), :isolated=>mean=>:isofraction)

# ╔═╡ ebe2194d-e79d-4a5e-bd83-a06eae95932f
isotab.unit

# ╔═╡ 9a2edeed-d33b-45e3-8dfc-1d03b967accc
F = @formula isofraction ~ 1 * Not(isofraction)

# ╔═╡ 3989392b-2f71-4613-ad85-d2fb828379d7


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

# ╔═╡ 301d385b-7879-445d-a903-772c10862750
widths

# ╔═╡ 6bdcf863-9946-4ca3-ab02-fa6aebe4b91d
# ╠═╡ show_logs = false
begin
    boundary = prop_str == "y-x" ? nothing : Dict("currentPathLength"=>(0,100))
    @time G = adaptive.get_grid(beh, props; widths, thresh, maxrad, radiusinc, boundary);
end

# ╔═╡ 9e635078-bfdb-41bf-8730-e08a968d5e71
md"""
Implied linear width of `maxrad[$(nanmaximum(G.radii))]` => $(nanmaximum(G.radii) * sqrt(2)) 
"""

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

# ╔═╡ 890fe951-19bb-4c9a-a905-d798bb36c57e
O = @time adaptive.get_occupancy(beh, G);

# ╔═╡ 592d79b4-edf6-4a0c-af73-1d2805d6410e
plot(O, clim=(0,0.01), ylims=ylim)

# ╔═╡ d7175827-7528-4cfe-bf3f-d9971f682f49
# ╠═╡ disabled = true
#=╠═╡
O
  ╠═╡ =#

# ╔═╡ bef016cd-26d1-4de8-a970-182fe2b92e88
# Test field abilities
@time multiunit = @time adaptive.get_adaptivefield(spikes, G, O);
# @benchmark adaptive.get_adaptivefield(spikes, G, O);

# ╔═╡ 410128bf-2332-4aa2-91cb-441e0235cc4a
plot(multiunit)

# ╔═╡ b88c0ec1-b150-49be-828f-6c32bb770c48
begin
	@time units = adaptive.yartsev(spikes, G, O; widths=width, thresh, 
	                               filters=filts[:all]);
end;

# ╔═╡ 62a6b931-b4ef-4431-8e7f-14db6e011d00
shifted = Timeshift.shifted_fields(beh, spikes, shifts, G.props;
                               shiftbeh=false,
                               widths, 
							   filters=filts[datacut], 
							   thresh);

# ╔═╡ d0ed9a46-00bb-4ce2-a2db-50bc060ec976
SFs = Timeshift.ShiftedFields(shifted);

# ╔═╡ 7cb361f2-09c9-48ad-ad01-971ac273ecd3
begin
	f = types.matrixform(SFs);
	push_shiftmetric!(f, best_tau!; metric=:bitsperspike);
end

# ╔═╡ f4466cd6-9f84-4fe4-a420-e77efc0c8ff8
begin
	idx_typ = sortperm(f[:bestshift_bitsperspike][:,1])
	f_s = f[idx_typ,:]
	M=Matrix(f_s[:bitsperspike])
end

# ╔═╡ 6c9a2819-b121-4b9a-a94d-11b3eb85fa7c
size(M)

# ╔═╡ bc93bf68-b7ac-4926-a7fd-dfff5b657f6a
heatmap(hcat([Utils.norm_extrema(m) for m in  eachrow(M)]...)')

# ╔═╡ 79b0a225-9019-4379-b9d0-61ff3408a690
metrics.push_dims!(f_s)

# ╔═╡ 60ba0b7b-e30d-42c8-9ff6-e17b5163ff53
isotau = leftjoin( DataFrame([f_s[:bestshift_bitsperspike][:,1] f_s[:unit][:,1]], [:bestshift_bitsperspike, :unit]), isotab; on=:unit)

# ╔═╡ e1136496-d5fa-442b-85b0-a69379c0335a
@df isotau scatter(:bestshift_bitsperspike, :isofraction)

# ╔═╡ 3c7440bf-0e77-4287-870c-3b88b0539607
begin
	shifted_iso = @memoize Timeshift.shifted_fields(beh, @subset(spikes,:isolated .== true), shifts, G.props;
							   shiftbeh=false,
							   widths, 
							   filters=filts[datacut], 
							   thresh);
	
	shifted_adj = @memoize Timeshift.shifted_fields(beh, @subset(spikes,:isolated .== false), shifts, G.props;
							   shiftbeh=false,
							   widths, 
							   filters=filts[datacut], 
							   thresh);
end

# ╔═╡ 129a57c7-aa83-41dc-8071-d6d7213d56b9
begin
		is_sfs = ShiftedFields(shifted_iso);
		f_is = types.matrixform(is_sfs);
		push_shiftmetric!(f_is, best_tau!; metric=:bitsperspike);
end

# ╔═╡ a279c651-b6f2-46d1-a0d6-c529a64517bf
adj_sfs = ShiftedFields(shifted_adj);

# ╔═╡ 65bed4d7-d800-4664-a152-6408a3621ab0
f_adj = types.matrixform(adj_sfs);

# ╔═╡ 45ad932d-922f-4e54-862e-b5db0c2f3658
push_shiftmetric!(f_adj, best_tau!; metric=:bitsperspike);

# ╔═╡ f4da9ee3-9882-4cd3-b3df-7f67e5f60ce1
begin
	idx = sortperm(f_is[:bestshift_bitsperspike][:,1])
	f_sort = f_is[idx,:]
	M_iso=Matrix(f_sort[:bitsperspike])
	
	idx = sortperm(f_adj[:bestshift_bitsperspike][:,1])
	f_sort_adj = f_adj[idx,:]
	M_adj = Matrix(f_sort_adj[:bitsperspike])
end

# ╔═╡ 444a4624-c446-4995-b4c2-c46241d3c668
begin
	pppp=plot()
	labs = Dict(1=>"total", 2=>"isolated", 3=>"adjacent")
	for (i, x) in enumerate([f_s[:,1][:bestshift_bitsperspike], vec(f_is[:,1][:bestshift_bitsperspike]), vec(f_adj[:,1][:bestshift_bitsperspike])])
		violin!(fill(i, size(x)), x, c=:gray, label="")
		scatter!(fill(i, size(x)) .+ 0.1 .* randn(size(x)), x, label=labs[i])
	end
	pppp
end

# ╔═╡ aff279a1-e837-484c-a5bd-112331cf2769
begin
	pppp
	Plot.setfolder("nonlocality","isolated-shifted")
	Plot.save((;desc="shifted sum of isolated firing (percent normed, around median)", filt=datacut, norm=norm));
end

# ╔═╡ 1bb9e911-b027-4f2b-a6d1-bdf81d0f174b
begin
	ppppp=plot()
	for (i, x) in enumerate([vec(f_s[:,1][:bestshift_bitsperspike]), vec(f_is[:,1][:bestshift_bitsperspike]), vec(f_adj[:,1][:bestshift_bitsperspike])])
		plot!(ecdf(x), label=labs[i])
	end
	ppppp
end

# ╔═╡ b873baf8-ceaa-47b5-b805-76f6bcfa49b2
begin
	ppppp
	Plot.setfolder("nonlocality","isolated-shifted")
	Plot.save((;desc="ecdf compare shifted sum of isolated firing (percent normed, around median)", filt=datacut, norm=norm));
end

# ╔═╡ 18abe083-65fd-4b6f-9bee-589c57f16b3f
begin
	sh_iso,un_iso = getshifts(is_sfs),[x[1] for x in getunits(is_sfs)]
	sh_adj,un_adj = getshifts(adj_sfs),[x[1] for x in getunits(adj_sfs)]
end

# ╔═╡ 697ffe5f-d3e6-44a5-b8ec-57a04a84a4e4
begin

	h_typ = heatmap(shifts,1:length(units),hcat([Utils.norm_percent(m,.5) for m in eachrow(M)]...)'; clim=(-50,50))
	vline!([0], c=:black, linestyle=:dash, linewidth=2, legend=:none)
	
	h_iso = heatmap(sh_iso,1:length(un_iso),hcat([Utils.norm_percent(m,.5) for m in eachrow(M_iso)]...)'; clim=(-50,50))
	vline!([0], c=:black, linestyle=:dash, linewidth=2, legend=:none)
	
	h_adj = heatmap(sh_iso,1:length(un_adj),hcat([Utils.norm_percent(m,.5) for m in eachrow(M_adj)]...)'; clim=(-50,50))
	vline!([0], c=:black, linestyle=:dash, linewidth=2, legend=:none)

	L = Plots.@layout [[ddd{0.25w} orig_h dd{0.25w}]; grid(1,2)]
	b=Plot.create_blank_plot()
	hhhh=plot(b,h_typ,b, h_adj,h_iso; layout=L)
	
end

# ╔═╡ 1167b4a7-1185-46f3-8e78-0a2bd97903ae
begin
	hhhh
	Plot.setfolder("nonlocality","isolated-shifted")
	Plot.save((;desc="shifted heatmap of isolated firing (percent normed, around median)", filt=datacut, norm=norm));
end

# ╔═╡ 671cb2a4-e6c0-4c23-9e15-41c861d8f383
begin
	widths["stopWell"] = 0.50
	Utils.filtreg.register(beh,spikes,on="time",transfer=["stopWell"])
	shifted_wg = Timeshift.shifted_fields(beh, spikes, shifts, ["x","y","stopWell"];
	                               shiftbeh=false,
	                               widths, 
								   adaptive=false,
                                   metricfuncs=[metrics.bitsperspike,metrics.totalcount,metrics.maxrate,metrics.meanrate],
								   filters=filts[datacut], 
								   thresh);
end

# ╔═╡ bb23f597-480e-43d1-a527-71f05a426d1d
begin
	bps = OrderedDict(k=>v[:bitsperspike] for (k,v) in units)
	sort!(bps, by=k->bps[k], rev=true)
	units_ordered = [x[1] for x in keys(bps)]
end

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

# ╔═╡ 4d814c3e-97e1-491a-b1d8-c7ca9c628afd
μ_firing = begin
    Q = units[(;unit=unit)]
    nansum(reshape(Q.occ.prob, size(Q.occ.count)) .* Q.rate)
end

# ╔═╡ eefa56cc-f303-40ee-aa44-dc758eac750b
field = units[(;unit=unit)]

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

# ╔═╡ 16b22203-b96b-4899-80f7-6928864f0543
plot(SFs[unit,timeshift])

# ╔═╡ 39562c46-79fb-4a0c-b59e-099108258190
keys(filts)

# ╔═╡ 8fc086f1-e619-4d44-a22b-7e9b9feee7fc
keys(f[1])

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

# ╔═╡ 719c49c0-f156-472c-ba15-26719753621f


# ╔═╡ 42d6cb62-93a8-426f-9542-d66fb8dd4d80
using Memoization

# ╔═╡ f5f1dbb9-3623-4628-8407-8809cd3fb118
using Memoization

# ╔═╡ Cell order:
# ╟─ff1db172-c3ab-41ea-920c-1dbf831c1336
# ╟─0be7ba01-a316-41ce-8df3-a5ae028c74e7
# ╟─37d7f4fd-80a7-47d0-8912-7f002620109f
# ╠═44dde9e4-f9ca-11ec-1348-d968780f671c
# ╠═da7809bb-a94f-440e-96c1-02e1feae9fc3
# ╠═42d6cb62-93a8-426f-9542-d66fb8dd4d80
# ╠═5f6e31d3-7101-49aa-a289-39e3967aa3a8
# ╠═a1ad6173-5ead-4559-bddb-9aee6119b9d0
# ╟─2f8ac703-417c-4360-a619-e799d8bb594f
# ╟─31082fe7-ed61-4d37-a025-77420da3f24a
# ╟─d51ce0f2-03bf-4c88-9302-9fd4bc8621eb
# ╟─ff355ad4-42da-4493-ae56-3bc9f0d8627c
# ╠═cbfbe9c9-f7c5-4219-8d81-b335fe5f5ed6
# ╟─efcdc2f1-5e26-4534-953e-defae4bd8603
# ╟─301d385b-7879-445d-a903-772c10862750
# ╟─6bdcf863-9946-4ca3-ab02-fa6aebe4b91d
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
# ╠═f5f1dbb9-3623-4628-8407-8809cd3fb118
# ╠═c34cfe1f-3663-4901-a258-0015cfa74a38
# ╠═a0d00fd1-c587-4bee-aa14-d0366a7f65ae
# ╟─7c6cfeb1-2c78-4480-852b-aa06cc818f76
# ╠═44abcbd4-5f71-4924-b77d-9680cc96044f
# ╠═3f6cec15-f43d-43d2-866b-45a04ea5715f
# ╠═f2fc5e48-3adc-41d0-b3e9-36bca362415b
# ╟─4d814c3e-97e1-491a-b1d8-c7ca9c628afd
# ╠═eefa56cc-f303-40ee-aa44-dc758eac750b
# ╟─f9378d49-2f86-4088-bc6d-3b5b227b7c66
# ╠═d1fe3b02-34f5-46ae-8da8-4bac71c86d84
# ╠═62a6b931-b4ef-4431-8e7f-14db6e011d00
# ╠═014d81aa-5441-41c1-8b71-a347dd03ac9a
# ╠═d0ed9a46-00bb-4ce2-a2db-50bc060ec976
# ╠═316dab79-015f-46e4-88f5-7051529484e5
# ╠═16b22203-b96b-4899-80f7-6928864f0543
# ╠═7cb361f2-09c9-48ad-ad01-971ac273ecd3
# ╠═f4466cd6-9f84-4fe4-a420-e77efc0c8ff8
# ╠═6c9a2819-b121-4b9a-a94d-11b3eb85fa7c
# ╠═bc93bf68-b7ac-4926-a7fd-dfff5b657f6a
# ╠═69d37df3-69d3-4c91-8adf-ff8ecd4c9df1
# ╠═664fba69-e498-462a-8e84-522a22f2498d
# ╠═6390cdc5-7d0e-456b-a46b-359ef1bdc63d
# ╠═bf12b73f-c146-4ab2-be79-5b8097add3f3
# ╠═00936498-8061-43da-bdb3-0c7b000770bb
# ╠═ffe4b420-41ae-4c58-8bd6-d5d684f48b78
# ╠═fff27c22-7c2c-45d0-a670-6920b9bd6a31
# ╠═28b0690b-e491-4cb0-be43-a3d23fc4903a
# ╠═39562c46-79fb-4a0c-b59e-099108258190
# ╠═3c7440bf-0e77-4287-870c-3b88b0539607
# ╠═129a57c7-aa83-41dc-8071-d6d7213d56b9
# ╠═a279c651-b6f2-46d1-a0d6-c529a64517bf
# ╠═65bed4d7-d800-4664-a152-6408a3621ab0
# ╠═45ad932d-922f-4e54-862e-b5db0c2f3658
# ╠═8fc086f1-e619-4d44-a22b-7e9b9feee7fc
# ╠═f4da9ee3-9882-4cd3-b3df-7f67e5f60ce1
# ╠═18abe083-65fd-4b6f-9bee-589c57f16b3f
# ╠═697ffe5f-d3e6-44a5-b8ec-57a04a84a4e4
# ╠═1167b4a7-1185-46f3-8e78-0a2bd97903ae
# ╠═444a4624-c446-4995-b4c2-c46241d3c668
# ╠═aff279a1-e837-484c-a5bd-112331cf2769
# ╠═f4ddcd9b-a940-40b3-8a31-25c1b2e103d9
# ╠═1bb9e911-b027-4f2b-a6d1-bdf81d0f174b
# ╠═b873baf8-ceaa-47b5-b805-76f6bcfa49b2
# ╠═85cf402a-eb7b-4d3c-972b-810ac1708632
# ╠═ebe2194d-e79d-4a5e-bd83-a06eae95932f
# ╠═79b0a225-9019-4379-b9d0-61ff3408a690
# ╠═1d93735a-1573-44cf-a650-dc639566a027
# ╠═60ba0b7b-e30d-42c8-9ff6-e17b5163ff53
# ╠═bb7c51b0-bac8-4aba-83df-30f301be4e65
# ╠═e1136496-d5fa-442b-85b0-a69379c0335a
# ╠═9bdea77f-22cc-414a-b675-3bd8bbe7cdf8
# ╠═9a2edeed-d33b-45e3-8dfc-1d03b967accc
# ╠═671cb2a4-e6c0-4c23-9e15-41c861d8f383
# ╠═3989392b-2f71-4613-ad85-d2fb828379d7
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
# ╟─97f2daad-1190-46a2-8c1a-288e0177a29b
# ╟─d48e9bc9-4d67-4e7b-a0f7-25b975813ccd
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
# ╠═719c49c0-f156-472c-ba15-26719753621f
