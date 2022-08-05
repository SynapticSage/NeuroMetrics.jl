### A Pluto.jl notebook ###
# v0.19.9

#> [frontmatter]
#> title = "Adaptive Field Shifting"
#> description = "Exploring field shift data with adaptive sampling routines"

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
	using DataFramesMeta
      using Statistics
      GC.gc()

	  using GoalFetchAnalysis
	  import Utils
	  import Timeshift
	  import Plot
	  using Field.metrics
	  using Timeshift.shiftmetrics	  
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

# ╔═╡ e11aae39-41ff-433f-8225-c5230c7c5e2e
using StatsPlots

# ╔═╡ ff1db172-c3ab-41ea-920c-1dbf831c1336
md"""
#### 🚀 **Shifted adaptive receptive shifts**

Purpose: Having established adaptive sampling works, and suspecting its answers may differ from fixed methods, we turn now to how those fields evolve over shifts. This notebook turns an eye to the field structure and metrics that quantify field structure over the shifts. 

**TODO list**
* metrics(shifted fields)
* shifted field plot recipes
"""

# ╔═╡ 0be7ba01-a316-41ce-8df3-a5ae028c74e7
PlutoUI.TableOfContents(title="Shifted Adaptive RFs" )

# ╔═╡ 37d7f4fd-80a7-47d0-8912-7f002620109f
md"""
# Preamble
Hidden cell below imports pacakges and sets up the environment.
"""

# ╔═╡ cb4d5494-24a6-4dfc-980b-23ec48fca7cc
F = Timeshift.load_fields();

# ╔═╡ 42fe7ad7-68c6-43ac-9588-17b902b6891b
import Utils: filtreg

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

# ╔═╡ b2436c80-7290-416f-87c9-137cfb601588
md"""
# 🔑 Key selection 👇
First, we either select a previous key or run the calculation"""

# ╔═╡ cbddd54f-83a6-433a-8e0f-67116b476e2e
md"""
*New method of picking keys*

In order to make radio buttons to pick our keys, we need to know the structure of the key options...
"""

# ╔═╡ c20dd185-4a34-433d-9775-f88475514add
begin

	# Determine which properties CAN be toggled
	allsets   = collect(keys(F))
	totalkeys = union(keys.(allsets)...)
	uvals = OrderedDict()
	for key in totalkeys
		push!(uvals, key=>unique([getindex(s, key) for s in allsets if key ∈ propertynames(s)]))
	end

	# Determine how they would be represented as controls and actual objects
	actuals  = OrderedDict()
	controls = OrderedDict()
	rep(x) = replace(x, "\""=>"", "]"=>"", "["=>"")
	for (i,(K,V)) in enumerate(uvals)
		sortV = [sort(V)..., nothing]
		push!(controls, K => [rep("$v") for v in sortV])
		for v in sortV
			push!(actuals,  ((v === nothing) ? "nothing" : "$v") => v)
		end
	end
end

# ╔═╡ 9f733213-1c8d-44e1-9f65-2c9861801a29
md"These encode the possible control knob strings"

# ╔═╡ 14daafe3-9df8-491d-a7ec-004fa7761c6b
controls

# ╔═╡ b3ef91fd-dfd0-44e6-b5f3-c3ff9b61e4e6
md"And the value that the actual knobs would point to"

# ╔═╡ f78cb175-ebc5-43e9-aa3b-50e7b6a8b95d
actuals

# ╔═╡ a5932740-2295-4973-8d21-31ad8048bea2
md"And here is the actual control console created from that control datastructure"

# ╔═╡ c461e634-d011-49bb-9ba5-b01df547f36f
begin
	datacut_sel = @bind datacut_key PlutoUI.Radio(controls[:datacut], default="all")
	thresh_sel  = @bind thresh_key PlutoUI.Radio(controls[:thresh], default="1.5")
	widths_sel  = @bind widths_key PlutoUI.Radio(controls[:widths], default="5.0")
	coact_sel   = @bind coactive PlutoUI.Radio(controls[:coactivity], default="nothing")
	console = (;datacut_sel, thresh_sel, widths_sel, coact_sel)
end

# ╔═╡ 36874cc8-8036-4d9b-acd4-a1ee34f238de
md"User selection is used to find the closest matching key in the cache"

# ╔═╡ d423498d-d024-4465-8dae-9eb899a75457
begin
	hard_settings = (;grid=:adaptive, first=-2.0, step=0.05, last=2.0)
	find_this = (;hard_settings..., coactivity=actuals[coactive], datacut=actuals[datacut_key], 			
				 thresh=actuals[thresh_key], widths=actuals[widths_key])
	key = Utils.namedtup.bestpartialmatch(keys(F), find_this;
										  nothing_means_removekey=true)
	#)
end

# ╔═╡ 0df2c833-c5f5-4664-b51a-8fc819d5a5f7
using Serialization; serialize(datadir("key.serial"), 
							   (;find_this, K=collect(keys(F)))
)

# ╔═╡ afcb55f9-ed2d-4860-a997-c272bece208f
md"grab that key, or if no key, compute with our settings above"

# ╔═╡ 955c7b75-00d4-4116-9242-92a7df8a0f87
shifted = F[key];

# ╔═╡ 832d1811-e3f8-492b-b44a-cc04edbd3a2c
begin
	pop_folder = Plot.setfolder("timeshift","population")
	md"""
	# 🫂 Population
	In this section, we take the complete set of shifts, and try to understand what they mean as a group.

	Data here is stored at $pop_folder
	"""
end

# ╔═╡ 2331218a-bc76-4adf-82e3-8e5b52aef0ca
plot_obj = Timeshift.DictOfShiftOfUnit{Float64}(shifted);

# ╔═╡ 5f15dc20-cf30-4088-a173-9c084ac2809a
# ╠═╡ show_logs = false
@time SFs = Timeshift.ShiftedFields(plot_obj);

# ╔═╡ 7946e555-1c73-4d13-a8f8-71ad6a2559f0
serialize(datadir("shiftedfields.serial"), SFs)

# ╔═╡ 0c392338-976c-4a5e-9bef-b15beb50ea22
cellfilters = begin
	kws=(;show_value=true)
	sl_spikecount = @bind spikecount Slider(50:10:150; default=50, kws...)
	sl_coherence  = @bind coherence Slider(0:0.1:1; default=0.6, kws...)
	sl_bitsperspike_ca1 = @bind bitsperspike_ca1 Slider(0:0.1:3; default=0.5, kws...)
	sl_bitsperspike_pfc = @bind bitsperspike_pfc Slider(0:0.05:3; default=0.25, kws...)
	(;sl_spikecount, sl_coherence, sl_bitsperspike_ca1, sl_bitsperspike_pfc)
end

# ╔═╡ b6208af9-4b51-47da-a01e-99773e87b853
begin
	(spikecount,) # doing this to signal dependency graph
_, prefilter = filtreg.register(cells, copy(SFs.metrics), on="unit", transfer=["area"]);
	prefilter
end

# ╔═╡ 3936c705-f0f1-45bb-902d-1c6bdb11ecc2
md"""
## 🤢 Filter out bad samples 🤮
In order for a neuron to survive to subsequent sections, it must survive these sieves. 
- Significance (from shuffle, if exists)
- Coherence (> current=$coherence), typical 0.7
- Information (> current=$bitsperspike_ca1), typical 0.5
- PFC spatial info (> current=$bitsperspike_pfc), lower than CA1, because these cells are not necessarily spatial
- Spikecount (> current=$spikecount), theoretical paper suggests 50 minimally
"""

# ╔═╡ bdee2a58-79e4-4053-a034-1db011687e66
begin
	sfsmet = copy(prefilter)
	filtrations(x) = ( # at least 1 shift passes this!
		(x.area .== "CA1" .&&
	     x.totalcount   .> spikecount .&&
	     x.coherence    .> coherence .&&
	     x.bitsperspike .> bitsperspike_ca1) .||
		(x.area .== "PFC" .&&
		 x.totalcount .> spikecount .&&
		 x.bitsperspike .> bitsperspike_pfc)
	)
	Filt.groupby_summary_condition_column!(sfsmet, :unit, filtrations, :area,
	                                       :totalcount, :coherence, :bitsperspike)
	deleteat!(sfsmet, sfsmet.condition .!= true)
end;

# ╔═╡ 166bffb7-72b9-4e9c-aa64-33db95763310
md"The full blow list of passing cells"

# ╔═╡ d9c1e7c8-f0b4-4bf6-afb4-27144ffa39f0
passingcells = combine(groupby(sfsmet, :unit), :area => first => :area, [:totalcount, :coherence, :bitsperspike] .=> maximum)

# ╔═╡ 1b36008f-4619-4fec-863c-abed47676e27
md"""
So, given the criteria our current data split (split=$datacut_key), we accept the following count of cells
"""

# ╔═╡ a9422422-65e3-415e-8404-acf9d8ed0e8a
begin
cs, ps = combine(groupby(cells, :area),  nrow),
combine(groupby(passingcells, :area),  nrow)
(;before=cs,after=ps)
end

# ╔═╡ f00d7750-ea75-403a-b221-49c8f338cec8
cellfilters

# ╔═╡ 2ad92fcd-48aa-43c7-b290-f4a5755592fd
begin
before = bar(cs.area, cs.nrow, c=:black,label="total")
after =  bar(ps.area, ps.nrow, c=:red,label="passing cells", title="cell filter results\n key=$datacut_key")
pfc_accepted = ps[ps.area .== "PFC",:nrow][1]
ca1_accepted = ps[ps.area .== "CA1",:nrow][1]
annotate!("CA1", ca1_accepted+5, text(ca1_accepted, :white))
annotate!("PFC", pfc_accepted+5, text(pfc_accepted, :white))
plot(before, after, link=:y)
end

# ╔═╡ 09a52111-e90f-4c10-a9ac-290abc489d2b
md"""
## Add sortable properties
"""

# ╔═╡ 8e6e4436-ceaf-432e-9a26-59de0642df96
for (on, measure) in Iterators.product([:bitsperspike, :coherence, :maxrate],			 								   [:best_tau!, :worst_tau!])
	shiftmetrics.metricapply!(sfsmet, measure; metric=on)
end;

# ╔═╡ 02c1873c-b59a-40a0-846d-1408b170cf37
md"Before plotting, this is the table plots will draw from"

# ╔═╡ 097466a4-1fe2-4ca5-bac3-fbc9ae53bc7e
sfsmet

# ╔═╡ 2cb7b7e0-e8aa-469c-833b-1a63ab21981a
begin
	sort_met_func(k) = occursin("best", string(k)) || occursin("worst", string(k))
	
	metric_sel = @bind met PlutoUI.Radio([String(k)=>k for k in propertynames(sfsmet)
						if k ∉ [:unit,:shift] && !(sort_met_func(k))],
				  default="bitsperspike")
	sort_sel   = @bind srt PlutoUI.Radio([String(k)=>k for k in 			
									propertynames(sfsmet)
									if sort_met_func(k)],
				 					default="bestshift_bitsperspike")
	
	desc_heat = @bind desc_heatmap PlutoUI.TextField();
	norm_by_range = @bind normrange PlutoUI.CheckBox()
	save_on = @bind allow_save_heat PlutoUI.CheckBox(default=true)
	(;datacut_sel, metric_sel, sort_sel, norm_by_range, saving=(;desc_heat, save_on))
end

# ╔═╡ 4dfc9b0e-234a-4f02-bbe8-1d3ff0f8dc42
md"""
## Cell sorted individual metrics
current: _$met_ by _$srt_
"""

# ╔═╡ 85d0cdd7-1154-4ec5-a35d-91080e54c415
md"""
Here, we calculate a stat matrix of **$met** sorted by **$srt**
"""

# ╔═╡ 90b33a8a-87fd-4d56-a609-41477dc8d93d
begin
	tmpstatmat = shiftmetrics.getstatmat(sfsmet, met; filtval=(met == "coherence" ? NaN : 0), asmat=true, unitnoskip=false, sortby=[srt], (normrange ? (;rangenorm=[0,1]) : (;percentnorm=0.5))...);
	
	plot_heatdist = heatmap(tmpstatmat.axes[2].val, 1:length(tmpstatmat.axes[1].val), tmpstatmat, ylabel="unit", title="$met by $srt\n$desc_heatmap\n", clim=(normrange ? (-Inf,Inf) : (-25,50)), colorbar_title=(normrange ? "Neuron Min(0) to Max(1)" : "Percent above median per neuron"), colorbar_titlefontrotation= 180)
	vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")

	annotate_pfc = true
	if annotate_pfc
		sfsmet_pfc = sort(@subset(sfsmet, :area .== "PFC", :shift .== 0), srt)
		rows_to_mark = findall(any(sfsmet_pfc.unit' .∈ tmpstatmat.axes[1].val; dims=2)[:,1])
		annotate!(sfsmet_pfc[:,srt], rows_to_mark, text("*", :white))
	end
	
	shifts = Timeshift.types.getshifts(shifted)
	xlim = (minimum(shifts), maximum(shifts))
	
	plot_histdist = histogram(sfsmet[:,srt],  normalize=:pdf, alpha=0.5, label="", xlabel="shift", ylabel="pdf")
	density!(sfsmet[:,srt]; xlim, label="")
	vline!([0], c=:black, linestyle=:dash, linewidth=2, label="")

	blank = plot(legend=false,grid=false,framestyle=:none,background_color_inside=:match)

	lay = @layout [a;b c{0.015w}]
	plot_shift_pop_overall = plot(plot_heatdist, plot_histdist, blank, layout=lay, link=:x, size=(700,700))
	
end

# ╔═╡ 18739638-ccac-4f1f-ac32-34c7c48bfb9c
allow_save_heat ? md"Saving..." : md"`save_on = false`"

# ╔═╡ 2fa18882-7e03-45a3-b794-cb44fec685b2
# ╠═╡ show_logs = false
if allow_save_heat
	plot(plot_shift_pop_overall)
	norm_heat_name = normrange ? "range" : "percent"
	Plot.save((;datacut=datacut_key, metric=met, sort=srt, norm=norm_heat_name, ((desc_heatmap != "") ? (;desc=desc_heatmap) : (;))...); rmexist=(:desc));
end

# ╔═╡ 5a67c4ae-1fc8-46b4-bba9-4ce3844feac4
begin
	intermetric_folder = Plot.setfolder("timeshift","population","intermetric-structure")
	md"""
	## No inter-metric structure
	(data stored in $intermetric_folder)
	"""
end

# ╔═╡ c8418c36-bac7-4e42-ae84-23809eb11642
label_unit_select = @bind label_unit CheckBox(default=true)

# ╔═╡ 6a121b95-67f9-4457-997d-9ea6efc270c7
begin
sfs_im = copy(sfsmet)
colorscheme = label_unit ? :glasbey_bw_minc_20_n256 : :vik 
color_field = label_unit ? :unit : :shift
colors = get(colorschemes[colorscheme], 					
			 Utils.norm_extrema(sfsmet[:,color_field],[0,1]))
(;colorscheme, color_field)
end

# ╔═╡ 523e49aa-4d75-4d2b-ac80-2502766a0705
begin
	#delete_neuron
	button_delete_neuron = @bind delete_neuron LabelButton("Delete neuron")
	button_shuffle_neurons = @bind shuffle_neurons LabelButton("Shuffle neurons")
	(button_delete_neuron, button_shuffle_neurons)
end

# ╔═╡ c5acc389-4274-4e47-9bcf-6e61799d34cf
begin
	delete_neuron
	del_msg = md"deleting neuron = $(sfs_im.unit[1])"
	@subset!(sfs_im, :unit .!= sfs_im.unit[1])
	del_msg
end

# ╔═╡ f61f2d87-54f8-4e30-9bd8-b96021ecdf36
# ╠═╡ show_logs = false
md"""
Two common metrics for examining what *is* a place field are **spatial coherence** and **bits per spike**. It's perfectly natural to question if they're providing similar or different information about the field.  Scattering each cells' shift measurements, suddenly we're very clear on that. These are for the most part not correlated metrics. That said, they do seem to exist on a $y=1/x$ manifiold, so they're not totally unrelated.
"""

# ╔═╡ 9b81e58c-9974-4261-b2b1-400f4954bf03
# ╠═╡ show_logs = false
begin
	delete_neuron
	p_scat_bits_coh =@df sfs_im scatter(:bitsperspike, :coherence, c=colors, xlabel="bitsperspike",label="")
	xx = 0.01:0.02:30; yy=1 ./ xx;
	plot!(xx,yy,c=:white,linestyle=:dash, linewidth=2, label="y=1/x")
	ylims!(0,2.5)
	xlims!(0,30)
	p_scat_bits_coh
end

# ╔═╡ b0055d27-682c-42d6-b086-6486c6889786
md"""
With the $y=1/x$ in mind, an alternative metric could be how far points move **normal** to $y=/1x$ (mahalanobis). And that *would* measures *how one moves against this determinstic relation* in the direction of organized and high information fields. This could be established by creating a separate $y=1/x$ list of points, and picking each cells' farthest shift from this $y=1/x$ list. 

It's also natural to wonder if this is modulated by the max firing rate of the field and these two metrics. Indeed, coherence seems to be every so slightly correlated with coherence. But not so much bits per spike. Although it is true that highest bits per spike confine to low max firing rate.
"""

# ╔═╡ 7bcd8e0e-e4bb-458c-9591-c51aec5f5c44
# ╠═╡ show_logs = false
# ╠═╡ disabled = true
#=╠═╡
@gif for azim in 1:1:360
	@df sfsmet scatter(:bitsperspike, :coherence, :maxrate; c=get(ColorSchemes.glasbey_bw_minc_20_n256, Utils.norm_extrema(:unit,[0,1])), camera=(azim,30),
		xlabel="bitsperspike", ylabel="coherence", zlabel="maxrate", label=""
	)
end
  ╠═╡ =#

# ╔═╡ f4a2ea03-e4bb-482e-95be-e2ddf4b38662
md"""
With either metric, there is no relationship in the best τ. Spatial coherence wouldn't identify the same best shifts. This could be suspected from prior plots.
"""

# ╔═╡ ec23bc9b-845c-46e5-88b2-99d900fc808b
@df sfsmet scatter(:bestshift_bitsperspike, :bestshift_coherence, c=get(ColorSchemes.glasbey_bw_minc_20_n256, Utils.norm_extrema(:unit,[0,1])), xlabel="τ_bitsperspike", ylabel="τ_coherence",label="")

# ╔═╡ 2cb37570-42f0-4136-a124-a8b7e2f5af62
md"""
Nor do any of these shifts correspond to picking shifts with higher max rates.
"""

# ╔═╡ 9133177f-1a80-4121-8615-56d4044e7bf6
# ╠═╡ disabled = true
#=╠═╡
@gif for azim in 1:1:360
	@df sfsmet scatter(:bestshift_bitsperspike, :bestshift_coherence, :bestshift_maxrate; c=get(ColorSchemes.glasbey_bw_minc_20_n256, Utils.norm_extrema(:unit,[0,1])), camera=(azim,30),
		xlabel="τ_{bitsperspike}", ylabel="τ_{coherence}", zlabel="τ_{maxrate}",
		label=""
	)
end
  ╠═╡ =#

# ╔═╡ fa3ec9e8-35c8-407b-8939-4e9d557d448b
# ╠═╡ disabled = true
#=╠═╡
@df SFs.metrics[:, [:bitsperspike, :coherence, :maxrate]] cornerplot(cols(1:3), compact=true)
  ╠═╡ =#

# ╔═╡ be6048d8-6f30-4d48-a755-5537c3b0b104
begin
sing_folder = Plot.setfolder("timeshift","single_cell")
md"""
# 🧜 Single field
This section is all about looking at the individual RFs that comprise the population above. This can be used to spot trouble.
Data is stored in $sing_folder
## View RF
Let's view a receptive field at a shift.
"""
end

# ╔═╡ 1994b7d7-c152-4aa7-a088-3f272246d9ae
# ╠═╡ disabled = true
#=╠═╡
nonzero_units = sort(passingcells.unit)
  ╠═╡ =#

# ╔═╡ 5731b325-0d82-4141-9122-3b67650ca2a5
md"Select out a unit 🦠 and a shift 🕘"

# ╔═╡ 5acf0a77-9e40-4117-83fa-4a0791849265
#=╠═╡
begin
	unit_sel = @bind shift_unit PlutoUI.Slider(nonzero_units, show_value=true)
	 shift_sel = @bind shift_shift PlutoUI.Slider(sort(collect(keys(shifted))), show_value=true,default=0)
	unit_shift_sel = (;unit_sel, shift_sel)
end
  ╠═╡ =#

# ╔═╡ 2a213e93-eaca-416d-9771-788e115e4081
#=╠═╡
begin
grabstr="And grab that data for that field (🦠 = $(shift_unit), 🕘 = $(shift_shift))"
md"$grabstr"
end
  ╠═╡ =#

# ╔═╡ 3afa8aae-bf3e-4364-8ad9-76fd50ca5ac9
#=╠═╡
single_shifted_field = get(plot_obj, shift_shift, shift_unit);
  ╠═╡ =#

# ╔═╡ 15e359aa-6ce2-4479-87cc-4ec3c7a5dfa2
md"Which we then plot"

# ╔═╡ dc714977-349e-417f-a5fd-1bd3a8c1e38b
md"""
key = $key
"""

# ╔═╡ 94930aab-8bb0-4da0-b26b-35ddb3efde3b
# ╠═╡ show_logs = false
#=╠═╡
begin
    plot(single_shifted_field; 
		aspect_ratio, ylims=ylim, 
		title=string(single_shifted_field.metrics)
	)
end
  ╠═╡ =#

# ╔═╡ dd9ad22d-25ff-4f6f-bf16-726c546a3829
md"And this panel can be reused here to compare with another dataset"

# ╔═╡ 89e8dc9a-b7d9-4d1f-b604-8376b28bfca5
#=╠═╡
(;console,unit_shift_sel)
  ╠═╡ =#

# ╔═╡ 47af1633-99bd-4dc2-9d91-9073ec327f27
md"""
## Plot properties
"""

# ╔═╡ ac9cddcb-097c-42a8-bd59-19fced23bf5a
# ╠═╡ disabled = true
#=╠═╡
shmet = sort(SFs[shift_unit].metrics,:shift);
  ╠═╡ =#

# ╔═╡ e070f95c-1671-4fd7-85b4-20b773de915e
#=╠═╡
(;console,unit_shift_sel)
  ╠═╡ =#

# ╔═╡ 8e872986-f947-4c0e-87bb-6e9db34a2508
md"""
Add spatial coherence to this
"""

# ╔═╡ 4e394823-6460-4aa7-af24-825a3f379a58
#=╠═╡
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

  ╠═╡ =#

# ╔═╡ d0be8e9b-5f26-46b0-96a0-7ba140b3135f
#=╠═╡
@df shmet scatter(:coherence, :bitsperspike, c=get(ColorSchemes.acton, Utils.norm_extrema(:shift,[0,1])), xlabel="Spatial Coherence", ylabel="Bits per spike", label="Unit = $shift_unit", legend=:topleft)
  ╠═╡ =#

# ╔═╡ 123070a2-0815-49f5-855a-2b87a120e16b
#=╠═╡
@gif for azim in 1:360
	@df shmet scatter3d(:shift, :coherence, :bitsperspike, c=get(ColorSchemes.acton, Utils.norm_extrema(:shift,[0,1])), xlabel="shift",ylabel="coherence", zlabel="bitsperspike", camera=(azim,30))
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─ff1db172-c3ab-41ea-920c-1dbf831c1336
# ╠═0be7ba01-a316-41ce-8df3-a5ae028c74e7
# ╟─37d7f4fd-80a7-47d0-8912-7f002620109f
# ╟─44dde9e4-f9ca-11ec-1348-d968780f671c
# ╠═cb4d5494-24a6-4dfc-980b-23ec48fca7cc
# ╠═42fe7ad7-68c6-43ac-9588-17b902b6891b
# ╟─a1ad6173-5ead-4559-bddb-9aee6119b9d0
# ╟─2f8ac703-417c-4360-a619-e799d8bb594f
# ╟─31082fe7-ed61-4d37-a025-77420da3f24a
# ╟─d51ce0f2-03bf-4c88-9302-9fd4bc8621eb
# ╟─ff355ad4-42da-4493-ae56-3bc9f0d8627c
# ╟─04eec95d-b5cf-47a5-a880-3146088cab00
# ╟─b2436c80-7290-416f-87c9-137cfb601588
# ╠═cbd7498a-5075-4e67-a829-95f9936146db
# ╟─cbddd54f-83a6-433a-8e0f-67116b476e2e
# ╟─c20dd185-4a34-433d-9775-f88475514add
# ╟─9f733213-1c8d-44e1-9f65-2c9861801a29
# ╟─14daafe3-9df8-491d-a7ec-004fa7761c6b
# ╟─b3ef91fd-dfd0-44e6-b5f3-c3ff9b61e4e6
# ╟─f78cb175-ebc5-43e9-aa3b-50e7b6a8b95d
# ╟─a5932740-2295-4973-8d21-31ad8048bea2
# ╟─c461e634-d011-49bb-9ba5-b01df547f36f
# ╟─36874cc8-8036-4d9b-acd4-a1ee34f238de
# ╟─d423498d-d024-4465-8dae-9eb899a75457
# ╠═0df2c833-c5f5-4664-b51a-8fc819d5a5f7
# ╟─afcb55f9-ed2d-4860-a997-c272bece208f
# ╠═955c7b75-00d4-4116-9242-92a7df8a0f87
# ╟─832d1811-e3f8-492b-b44a-cc04edbd3a2c
# ╠═2331218a-bc76-4adf-82e3-8e5b52aef0ca
# ╠═5f15dc20-cf30-4088-a173-9c084ac2809a
# ╠═7946e555-1c73-4d13-a8f8-71ad6a2559f0
# ╟─b6208af9-4b51-47da-a01e-99773e87b853
# ╟─3936c705-f0f1-45bb-902d-1c6bdb11ecc2
# ╟─0c392338-976c-4a5e-9bef-b15beb50ea22
# ╟─bdee2a58-79e4-4053-a034-1db011687e66
# ╟─166bffb7-72b9-4e9c-aa64-33db95763310
# ╟─d9c1e7c8-f0b4-4bf6-afb4-27144ffa39f0
# ╟─1b36008f-4619-4fec-863c-abed47676e27
# ╟─a9422422-65e3-415e-8404-acf9d8ed0e8a
# ╟─f00d7750-ea75-403a-b221-49c8f338cec8
# ╟─2ad92fcd-48aa-43c7-b290-f4a5755592fd
# ╟─09a52111-e90f-4c10-a9ac-290abc489d2b
# ╠═8e6e4436-ceaf-432e-9a26-59de0642df96
# ╟─4dfc9b0e-234a-4f02-bbe8-1d3ff0f8dc42
# ╟─85d0cdd7-1154-4ec5-a35d-91080e54c415
# ╟─02c1873c-b59a-40a0-846d-1408b170cf37
# ╟─097466a4-1fe2-4ca5-bac3-fbc9ae53bc7e
# ╟─2cb7b7e0-e8aa-469c-833b-1a63ab21981a
# ╟─90b33a8a-87fd-4d56-a609-41477dc8d93d
# ╟─18739638-ccac-4f1f-ac32-34c7c48bfb9c
# ╟─2fa18882-7e03-45a3-b794-cb44fec685b2
# ╟─5a67c4ae-1fc8-46b4-bba9-4ce3844feac4
# ╟─c8418c36-bac7-4e42-ae84-23809eb11642
# ╟─6a121b95-67f9-4457-997d-9ea6efc270c7
# ╟─523e49aa-4d75-4d2b-ac80-2502766a0705
# ╟─c5acc389-4274-4e47-9bcf-6e61799d34cf
# ╟─f61f2d87-54f8-4e30-9bd8-b96021ecdf36
# ╟─9b81e58c-9974-4261-b2b1-400f4954bf03
# ╟─b0055d27-682c-42d6-b086-6486c6889786
# ╠═7bcd8e0e-e4bb-458c-9591-c51aec5f5c44
# ╟─f4a2ea03-e4bb-482e-95be-e2ddf4b38662
# ╠═ec23bc9b-845c-46e5-88b2-99d900fc808b
# ╟─2cb37570-42f0-4136-a124-a8b7e2f5af62
# ╠═9133177f-1a80-4121-8615-56d4044e7bf6
# ╟─fa3ec9e8-35c8-407b-8939-4e9d557d448b
# ╟─be6048d8-6f30-4d48-a755-5537c3b0b104
# ╠═1994b7d7-c152-4aa7-a088-3f272246d9ae
# ╟─5731b325-0d82-4141-9122-3b67650ca2a5
# ╟─5acf0a77-9e40-4117-83fa-4a0791849265
# ╟─2a213e93-eaca-416d-9771-788e115e4081
# ╠═3afa8aae-bf3e-4364-8ad9-76fd50ca5ac9
# ╟─15e359aa-6ce2-4479-87cc-4ec3c7a5dfa2
# ╟─dc714977-349e-417f-a5fd-1bd3a8c1e38b
# ╟─94930aab-8bb0-4da0-b26b-35ddb3efde3b
# ╟─dd9ad22d-25ff-4f6f-bf16-726c546a3829
# ╠═89e8dc9a-b7d9-4d1f-b604-8376b28bfca5
# ╟─47af1633-99bd-4dc2-9d91-9073ec327f27
# ╠═ac9cddcb-097c-42a8-bd59-19fced23bf5a
# ╠═e11aae39-41ff-433f-8225-c5230c7c5e2e
# ╟─e070f95c-1671-4fd7-85b4-20b773de915e
# ╟─8e872986-f947-4c0e-87bb-6e9db34a2508
# ╟─4e394823-6460-4aa7-af24-825a3f379a58
# ╟─d0be8e9b-5f26-46b0-96a0-7ba140b3135f
# ╟─123070a2-0815-49f5-855a-2b87a120e16b
