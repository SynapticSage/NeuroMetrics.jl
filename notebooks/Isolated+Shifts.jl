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

# ╔═╡ e39a9dbf-4249-404c-813a-38be1eb57c31
using Memoization

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

# ╔═╡ f9378d49-2f86-4088-bc6d-3b5b227b7c66
md"""
### 📜 `to_dataframe`

we're going to want codes that take a set of receptive fields and turn them into a dataframe

- fields
- metrics
"""

# ╔═╡ d1fe3b02-34f5-46ae-8da8-4bac71c86d84
shifts = -2:0.1:2

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
	spikeiso = true
end

# ╔═╡ bf12b73f-c146-4ab2-be79-5b8097add3f3
begin
	datacut_sl = @bind datacut_str PlutoUI.Radio(["all","cue","memory","nontask"], default="memory")

	norm_sl = @bind norm PlutoUI.Radio(["01","percent"],default="percent")
	
	iso_sl=(;datacut_sl, norm_sl)
end

# ╔═╡ 00936498-8061-43da-bdb3-0c7b000770bb
datacut = Symbol(datacut_str)

# ╔═╡ 62a6b931-b4ef-4431-8e7f-14db6e011d00
@time shifted = Timeshift.shifted_fields(beh, spikes, shifts, G.props;
                               shiftbeh=false,
                               widths, 
							   filters=filts[datacut], 
							   thresh);

# ╔═╡ d0ed9a46-00bb-4ce2-a2db-50bc060ec976
SFs = Timeshift.ShiftedFields(shifted);

# ╔═╡ 16b22203-b96b-4899-80f7-6928864f0543
begin
	Plots.@gif for sh in shifts
		plot(SFs[unit,sh])
	end every 1
end

# ╔═╡ 7cb361f2-09c9-48ad-ad01-971ac273ecd3
begin
	f = types.matrixform(SFs);
	push_metric!(f, metrics.bitsperspike)
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

# ╔═╡ ffe4b420-41ae-4c58-8bd6-d5d684f48b78
if norm == "percent"
		normf = Utils.norm_percent
	else
		normf = Utils.norm_extrema
	end

# ╔═╡ fff27c22-7c2c-45d0-a670-6920b9bd6a31
all((!).(spikes.isolated))

# ╔═╡ 39562c46-79fb-4a0c-b59e-099108258190
keys(filts)

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
		push_metric!(f_is, metrics.bitsperspike)
		push_shiftmetric!(f_is, best_tau!; metric=:bitsperspike);
end

# ╔═╡ a279c651-b6f2-46d1-a0d6-c529a64517bf
begin
	 adj_sfs = ShiftedFields(shifted_adj);
 	f_adj = types.matrixform(adj_sfs);
	push_metric!(f_adj, metrics.bitsperspike)
 	push_shiftmetric!(f_adj, best_tau!; metric=:bitsperspike);
end

# ╔═╡ 8fc086f1-e619-4d44-a22b-7e9b9feee7fc
keys(f[1])


# ╔═╡ f4da9ee3-9882-4cd3-b3df-7f67e5f60ce1
begin
	idx = sortperm(f_is[:bestshift_bitsperspike][:,1])
	f_sort = f_is[idx,:]
	M_iso=Matrix(f_sort[:bitsperspike])
	
	idx = sortperm(f_adj[:bestshift_bitsperspike][:,1])
	f_sort_adj = f_adj[idx,:]
	M_adj = Matrix(f_sort_adj[:bitsperspike])
end

# ╔═╡ 18abe083-65fd-4b6f-9bee-589c57f16b3f
begin
	sh_iso,un_iso = getshifts(is_sfs), [x[1] for x in getunits(is_sfs)]
	sh_adj,un_adj = getshifts(adj_sfs),[x[1] for x in getunits(adj_sfs)]
end

# ╔═╡ 5369cc71-1132-482c-9251-f61423f1e0f0
begin
	h_typ = heatmap(shifts,1:length(units),hcat([Utils.norm_percent(m,.5) for m in eachrow(M)]...)'; clim=(-50,50))
	vline!([0], c=:black, linestyle=:dash, linewidth=2, legend=:none)
	
	h_iso = heatmap(sh_iso,1:length(un_iso),hcat([Utils.norm_percent(m,.5) for m in eachrow(M_iso)]...)'; clim=(-50,50))
	vline!([0], c=:black, linestyle=:dash, linewidth=2, legend=:none)
	
	h_adj = heatmap(sh_iso,1:length(un_adj),hcat([Utils.norm_percent(m,.5) for m in eachrow(M_adj)]...)'; clim=(-50,50))
	vline!([0], c=:black, linestyle=:dash, linewidth=2, legend=:none)

	L = Plots.@layout [[ddd{0.25w} orig_h dd{0.25w}]; grid(1,2)]
	b=Plot.create_blank_plot()
	hhhh=plot(b,h_typ,b, h_adj,h_iso; layout=L, size=(500,1000))
	
end

# ╔═╡ 1167b4a7-1185-46f3-8e78-0a2bd97903ae
begin
	hhhh
	Plot.setfolder("nonlocality","isolated-shifted")
	Plot.save((;desc="shifted heatmap of isolated firing (percent normed, around median)", filt=datacut, norm=norm));
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

# ╔═╡ 85cf402a-eb7b-4d3c-972b-810ac1708632
isotab = combine(groupby(spikes, :unit), :isolated=>mean=>:isofraction)

# ╔═╡ ebe2194d-e79d-4a5e-bd83-a06eae95932f
isotab.unit

# ╔═╡ 79b0a225-9019-4379-b9d0-61ff3408a690
metrics.push_dims!(f_s)

# ╔═╡ 60ba0b7b-e30d-42c8-9ff6-e17b5163ff53
isotau = leftjoin( 
	DataFrame([f_s[:bestshift_bitsperspike][:,1] f_s[:unit][:,1]], [:bestshift_bitsperspike, :unit]), 
	isotab; 
	on=:unit)

# ╔═╡ e1136496-d5fa-442b-85b0-a69379c0335a
@df isotau scatter(:bestshift_bitsperspike, :isofraction)

# ╔═╡ 9a2edeed-d33b-45e3-8dfc-1d03b967accc
F = @formula isofraction ~ 1 * Not(isofraction)

# ╔═╡ 671cb2a4-e6c0-4c23-9e15-41c861d8f383
begin
	widths["stopWell"] = 0.5
	Utils.filtreg.register(beh,spikes,on="time",transfer=["stopWell"])
	shifted_wg = Timeshift.shifted_fields(beh, spikes, shifts, ["x","y","stopWell"];
	                               shiftbeh=false,
	                               widths, 
								   adpative=false,
								   filters=filts[datacut], 
								   thresh);
end;

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
# ╟─f9378d49-2f86-4088-bc6d-3b5b227b7c66
# ╠═d1fe3b02-34f5-46ae-8da8-4bac71c86d84
# ╠═62a6b931-b4ef-4431-8e7f-14db6e011d00
# ╠═d0ed9a46-00bb-4ce2-a2db-50bc060ec976
# ╠═316dab79-015f-46e4-88f5-7051529484e5
# ╠═16b22203-b96b-4899-80f7-6928864f0543
# ╠═7cb361f2-09c9-48ad-ad01-971ac273ecd3
# ╠═f4466cd6-9f84-4fe4-a420-e77efc0c8ff8
# ╠═6c9a2819-b121-4b9a-a94d-11b3eb85fa7c
# ╠═bc93bf68-b7ac-4926-a7fd-dfff5b657f6a
# ╠═69d37df3-69d3-4c91-8adf-ff8ecd4c9df1
# ╠═664fba69-e498-462a-8e84-522a22f2498d
# ╠═bf12b73f-c146-4ab2-be79-5b8097add3f3
# ╠═00936498-8061-43da-bdb3-0c7b000770bb
# ╠═ffe4b420-41ae-4c58-8bd6-d5d684f48b78
# ╠═fff27c22-7c2c-45d0-a670-6920b9bd6a31
# ╠═39562c46-79fb-4a0c-b59e-099108258190
# ╠═e39a9dbf-4249-404c-813a-38be1eb57c31
# ╠═3c7440bf-0e77-4287-870c-3b88b0539607
# ╠═129a57c7-aa83-41dc-8071-d6d7213d56b9
# ╠═a279c651-b6f2-46d1-a0d6-c529a64517bf
# ╠═8fc086f1-e619-4d44-a22b-7e9b9feee7fc
# ╠═f4da9ee3-9882-4cd3-b3df-7f67e5f60ce1
# ╠═18abe083-65fd-4b6f-9bee-589c57f16b3f
# ╠═5369cc71-1132-482c-9251-f61423f1e0f0
# ╠═1167b4a7-1185-46f3-8e78-0a2bd97903ae
# ╠═444a4624-c446-4995-b4c2-c46241d3c668
# ╠═aff279a1-e837-484c-a5bd-112331cf2769
# ╠═1bb9e911-b027-4f2b-a6d1-bdf81d0f174b
# ╠═b873baf8-ceaa-47b5-b805-76f6bcfa49b2
# ╠═85cf402a-eb7b-4d3c-972b-810ac1708632
# ╠═ebe2194d-e79d-4a5e-bd83-a06eae95932f
# ╠═79b0a225-9019-4379-b9d0-61ff3408a690
# ╠═60ba0b7b-e30d-42c8-9ff6-e17b5163ff53
# ╠═e1136496-d5fa-442b-85b0-a69379c0335a
# ╠═9a2edeed-d33b-45e3-8dfc-1d03b967accc
# ╠═671cb2a4-e6c0-4c23-9e15-41c861d8f383
