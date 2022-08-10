### A Pluto.jl notebook ###
# v0.19.9

#> [frontmatter]
#> title = "Timeshift: Compare mains"

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

# ╔═╡ e8c8870a-0271-11ed-2ac4-317a38722303
# ╠═╡ show_logs = false
begin
      using DrWatson
      quickactivate(expanduser("~/Projects/goal-code"))
      using Plots
      using Revise
      using DataFrames, DataFramesMeta
      using NaNStatistics
      import ProgressLogging
      using PlutoUI
      using DataStructures: OrderedDict
	  using Statistics, NaNStatistics
	  using StatsBase

      using GoalFetchAnalysis     
      import Utils 
      import Timeshift
      using Timeshift
      import Plot
	  using Plot.timeshift
	  using Utils.namedtup
	  using Utils: filtreg

      adaptive = Field.adaptive
      metrics = Field.metrics
      WIDTHS = OrderedDict(
          "x"=>2.5f0, "y"=>2.5f0, "currentPathLength"=>1f0,
          "currentAngle"=>Float32(2pi/80)
      )
end;

# ╔═╡ 9ce4497d-2b82-464b-9df9-c2e568a661ac
using Serialization

# ╔═╡ 51e1b371-da12-403a-ba4e-ce7790bd5d2a
using StatsPlots

# ╔═╡ 5005402b-4d25-41a3-916b-4c814faa9065
md"""
#### Compare Dataset Metrics
Purupose: This notebook exists to do some quick comparisons between datasets. Not intended to be an analysis workhorse.
"""

# ╔═╡ 44079721-eb63-4c82-af86-ee3b4f55ab42
PlutoUI.TableOfContents(title="Compare Dataset Metrics", aside=true)

# ╔═╡ 5c008e3c-9395-4f43-be60-ff510a54b97e
md"Loadup our libraries..."

# ╔═╡ e4174f16-1628-4e7c-8d30-840690422582
md"""
# Load
"""

# ╔═╡ e6605880-9bd0-43b3-81dd-643829d32cbc
cells = Load.load_cells("RY16",36);

# ╔═╡ 99b0eff5-6d22-4948-91ae-5080d838580e
@time I = Timeshift.load_mains();

# ╔═╡ ffa1cf80-e6c5-4d5f-b3fb-555018fc5f75
function table_transform(x)
	x = transform(Table.to_dataframe( x, key_name=[:shift], explode=true), All(), :dim_1 => (x->replace(x,missing=>0)), renamecols=false)
	x = flatten(x, [:value, :dim_1])
	x = unstack(x, :metric, :value)
	x = filtreg.register(cells, x, on="unit", transfer=["area"])[2]

	dropmissing!(x, :bitsperspike)
	for (on, measure) in Iterators.product([:bitsperspike, :coherence, :maxrate],			 								   [:best_tau!, :worst_tau!])
		shiftmetrics.metricapply!(x, measure; metric=on)
	end

	x
end

# ╔═╡ 01000490-e094-4c88-ac60-7b9fdeba3ccc
keys(I)

# ╔═╡ c753704b-f9bc-4824-a3ec-27b6f8ff3f8b
begin
    controls, actuals = Timeshift.pluto_keyselector(I)
    md"""
    ## Keyset 1
    """
end

# ╔═╡ 502701a1-9c58-48ca-b85b-b580b6fecde7
begin
	datacut_sel1 = @bind datacut_key1 PlutoUI.Radio(controls[:datacut], default="all")
	thresh_sel1  = @bind thresh_key1 PlutoUI.Radio(controls[:thresh], default="1.5")
	widths_sel1  = @bind widths_key1 PlutoUI.Radio(controls[:widths], default="5.0")
	#coact_sel1   = @bind coactive1 PlutoUI.Radio(controls[:coactivity], default="nothing")
	console1 = (;datacut_sel1, thresh_sel1, widths_sel1)#, coact_sel1)

end

# ╔═╡ 635d367d-accb-438d-b67d-3a9e5fdd7cb7
key1 = Timeshift.pluto_executekey(I, actuals, :thresh=>thresh_key1, 
                                  :widths=>widths_key1, 
                                  #:coactivity=>coactive1, 
                                  :datacut=>datacut_key1)

# ╔═╡ 63772ff9-2b0c-4d30-a084-3369e3e809e0
# ╠═╡ show_logs = false
begin
	sf1 = table_transform(I[key1])
	@assert length(unique(sf1.bestshift_bitsperspike)) > 1
	sf1
end;

# ╔═╡ 4cfc656c-957c-4394-833d-93af05eff81c
md"""
## Keyset 2
"""

# ╔═╡ d600f8bf-555b-4972-94b3-fc0a07a241c0
begin
	datacut_sel2 = @bind datacut_key2 PlutoUI.Radio(controls[:datacut], default="all")
	thresh_sel2  = @bind thresh_key2 PlutoUI.Radio(controls[:thresh], default="1.5")
	widths_sel2  = @bind widths_key2 PlutoUI.Radio(controls[:widths], default="5.0")
	#coact_sel2   = @bind coactive2 PlutoUI.Radio(controls[:coactivity], default="nothing")
	console2 = (;datacut_sel2, thresh_sel2, widths_sel2)#, coact_sel2)
end

# ╔═╡ 21e1d7d4-0120-41b6-81c0-f98819e0dd13
key2 = Timeshift.pluto_executekey(I, actuals, :thresh=>thresh_key2, 
                                  :widths=>widths_key2, 
                                  #:coactivity=>coactive2, 
                                  :datacut=>datacut_key2)

# ╔═╡ 4902ca06-9f29-48c7-a3f3-8a6178ddaa4c
# ╠═╡ show_logs = false
sf2 = table_transform(I[key2]);

# ╔═╡ 58591957-b16f-4c4a-b3e7-f6c36feca331
md"Shortcut names for keys"

# ╔═╡ c1053cca-23a9-40cc-a2db-bbbb07c8a928
begin
	f(K) = NamedTuple(filter(k-> k[1] ∉  [:first, :last, :step, :grid, :marginal], pairs(K)))
	k1, k2 = string(f(key1)), string(f(key2))
	k1f = replace(k1, ", "=>"\n")
	k2f = replace(k2, ", "=>"\n")
end;

# ╔═╡ 1a4cc678-6900-4e07-a6dc-78073df8b06a
k1

# ╔═╡ fb770592-cbed-4c47-b156-e2cdde09b529
begin
	sel_prop = @bind met Radio(["bitsperspike", "coherence", "maxrate"], default="bitsperspike")
	sel_sort = @bind srt Radio([string(x) for x in names(sf1) if occursin("best", x)], default="bestshift_bitsperspike")
	sel = (;prop=sel_prop, srt=sel_sort)
end

# ╔═╡ 41fbea69-ba6c-48be-8b81-244f2e01d428
md"""
## Filtration
filter datasets
"""

# ╔═╡ 9f47b399-bc14-4369-8432-7d2b8b5e872c
cellfilters = begin
	kws=(;show_value=true)
	sl_spikecount = @bind spikecount Slider(50:10:150; default=50, kws...)
	sl_coherence  = @bind coherence Slider(0:0.1:1; default=0.6, kws...)
	sl_bitsperspike_ca1 = @bind bitsperspike_ca1 Slider(0:0.1:3; default=0.5, kws...)
	sl_bitsperspike_pfc = @bind bitsperspike_pfc Slider(0:0.05:3; default=0.25, kws...)
	(;sl_spikecount, sl_coherence, sl_bitsperspike_ca1, sl_bitsperspike_pfc)
end

# ╔═╡ 4732f7ef-cec0-4435-a2aa-370babc06cd2
begin
	function filtercells(prefilter)
		sfsmet = copy(prefilter)
		dropmissing!(sfsmet, :bitsperspike)
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
	end
	sf1F, sf2F = filtercells(sf1), filtercells(sf2)
end;

# ╔═╡ ef5b86bb-57f9-464b-9cae-acbf8933de1d
passingcells(x) = combine(groupby(x, :unit), :area => first => :area, [:totalcount, :coherence, :bitsperspike] .=> maximum)

# ╔═╡ 2266bece-09f2-41a2-87dc-adfb16775361
passingcells(sf1F), passingcells(sf2F)

# ╔═╡ 404416fe-70d8-456a-90e7-52de79fa2eb5
md"""
# Examine
"""

# ╔═╡ 88fa5b09-0fab-4293-89a8-02ca6ba80aa0

md"""
## Metric Δdistribution
"""

# ╔═╡ 231d579e-23d6-4ac7-a96d-e37366277965
(;console1, console2)

# ╔═╡ 3079567e-7fb0-4269-a01e-74be600b459f
begin
	function metric_distribution(sf1,sf2)
		cs(x) = collect(Utils.skipnan(skipmissing(x)))
		q1 = quantile(cs(sf1[!,met]), 0.96)
		q2 = quantile(cs(sf2[!,met]), 0.96)
		xlim1=(0,q1)
		xlim2=(0,q2)
		xlim3=(0,max(q1,q2))
		e1 = ecdf(cs(sf1[!,met]))
		e2 = ecdf(cs(sf2[!,met]))
		pe=plot(e1, xlim=xlim1, label="$(key1.datacut)")
		plot!(e2, c=:yellow, xlim=xlim3, label="$(key2.datacut)", legend=:outerbottomright, )
	plot(
		pe,
		histogram(cs(sf1[!,met]), yscale=:log10, title=k1f, xlim=xlim3,
			xlabel=string(met), label="", linecolor=nothing),
		histogram(cs(sf2[!,met]), yscale=:log10, title=k2f, c=:yellow, xlim=xlim3,
			xlabel=string(met), label="", linecolor=nothing, alpha=0.5),
		link=:x,
		layout = @layout [a; b c]
)
	end
	metric_distribution(sf1,sf2)
end

# ╔═╡ 9fda323f-46e7-436d-a2a4-8094527b1b10
metric_distribution(sf1F,sf2F)

# ╔═╡ 5c962c04-b8a1-46e9-8933-b65c16d76ef1
md"""
## Plot
### Straightup side-by-side, unfitlered
"""

# ╔═╡ 4f98b564-6696-47f9-89cc-0c5b92137753
normrange=true

# ╔═╡ 81ff9d0e-6e15-4182-8dd9-ee10457348e7
md"Munge matrices encoding our snakes"

# ╔═╡ 8723bf8e-93cc-4d55-a274-7699daec4044
md"And encode the plot function"

# ╔═╡ 87171f23-3f8f-4c33-8430-2c5e13732675
begin
	function plot_sorted_mets(sfsmet, tmpstatmat; met, srt, normrange, key)

		plot_heatdist = heatmap(tmpstatmat.axes[2].val, 1:length(tmpstatmat.axes[1].val), tmpstatmat, ylabel="unit", title="$(key.datacut): $met by\n $srt\n", clim=(normrange ? (-Inf,Inf) : (-25,50)), colorbar_title=(normrange ? "Neuron Min(0) to Max(1)" : "Percent above median per neuron"), colorbar_titlefontrotation= 180)
		vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
		
		shifts = sort(unique(sfsmet.shift))
		#@info shifts
		xlim = (minimum(shifts), maximum(shifts))
		#@info xlim
		shiftdist = dropmissing(sfsmet, srt)[!,srt]
		#@info shiftdist

		annotate_pfc = true
		if annotate_pfc
			sfsmet_pfc = @subset(sfsmet, :area .== "PFC", :shift .== 0, :unit .∈ [tmpstatmat.axes[1].val])
			sort!(sfsmet_pfc, srt)
			rows_to_mark = findall(any(sfsmet_pfc.unit' .∈ tmpstatmat.axes[1].val; dims=2)[:,1])
			@assert issorted(sfsmet_pfc[:,srt])
			@assert issorted(rows_to_mark)
			@info sfsmet_pfc[:,srt]
			@assert size(rows_to_mark) == size(sfsmet_pfc[:,srt])
			annotate!(sfsmet_pfc[:,srt], rows_to_mark, text("*", :gray))
			serialize(datadir("test.serial"),(;sfsmet_pfc, rows_to_mark, tmpstatmat, srt))
		end
		
		plot_histdist = 
			histogram(shiftdist;  xlim, normalize=:pdf, alpha=0.5, label="", xlabel="shift", ylabel="pdf")
		vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
		density!(sfsmet[:,srt]; xlim, label="", ylim=(0,0.4))

		
		blank = plot(legend=false, grid=false, framestyle=:none, background_color_inside=:match)
	
		lay = @layout [a;b c{0.015w}]
		
		plot_shift_pop_overall = plot(plot_heatdist, plot_histdist, blank, layout=lay, link=:x, size=(700,700))
		
		plot_shift_pop_overall = plot(plot_heatdist, plot_histdist, blank, layout=lay, link=:x, size=(700,700))
		plots=plot_shift_pop_overall
	end
end

# ╔═╡ 1afbdcbf-3c8d-45ba-9101-8f1f0db13512
md"Which is run on each cut"

# ╔═╡ 7ce8f6db-6696-411e-8751-06dbcf1d25da
begin
	function get_metmatrix(sfsmet; met, srt, normrange)	
		tmpstatmat = shiftmetrics.getstatmat(sfsmet, met; filtval=(met == "coherence" ? NaN : 0), asmat=true, unitnoskip=false, sortby=[srt], (normrange ? (;rangenorm=[0,1]) : (;percentnorm=0.5))...);
	end
	function side_by_side_plot(sf1, sf2; met, srt, normrange)
		metmatrix1 = get_metmatrix(sf1; met, srt, normrange)
		metmatrix2 = get_metmatrix(sf2; met, srt, normrange)
		plot_shift_pop_overall1, plot_shift_pop_overall2 = plot_sorted_mets(sf1, metmatrix1; met, srt, normrange=true, key=key1), plot_sorted_mets(sf2, metmatrix2; met, srt, normrange=true, key=key2)
		plot(plot_shift_pop_overall1, plot_shift_pop_overall2, size=(800,400))
	end
	side_by_side_plot(sf1, sf2; met, srt, normrange)
end

# ╔═╡ 11f54afb-b588-49b9-81c1-8a8393084916
md"""
* Why are PFC asterisks in the wrong place?
"""

# ╔═╡ 53b84249-b032-4e76-9d35-6dd7856194d9
(;console1, console2)

# ╔═╡ b5f732e2-369f-4d34-b070-de4352592bec
md"""
### Straightup side-by-side, filtered
Filter the data
"""

# ╔═╡ b509332f-cd1d-4a6d-8b79-954b02be35b4
side_by_side_plot(sf1F, sf2F; met, srt, normrange)

# ╔═╡ a6b29b91-d0fc-46a3-8e12-47cfbb207184
key1.datacut

# ╔═╡ e24a55c4-2946-47a4-a12a-3ef6654865c9
key2.datacut

# ╔═╡ a0083bf5-a36c-4b70-a077-8758708721da
begin
	function futurepast(x)
		answer = Vector{Int}(undef, length(x))
		answer.= 0
		answer[x .> 0.33] .= 1
		answer[x .< -0.33] .= -1
		answer
	end
	transform!(sf1F, srt => futurepast => :futurepast)
	transform!(sf2F, srt => futurepast => :futurepast)
	sfF = vcat(sf1F, sf2F, cols=:union, source = :datacut => [string(key1.datacut), string(key2.datacut)])
end

# ╔═╡ f4baa229-883a-41df-ba42-581783029bb1
function futurepast_index(x; thresh=0.1)
	lower, upper = -thresh, thresh
	L, H = sum(x .< lower), sum(x .> upper)
	@info "run" L H lower upper x
	( (H - L) / (H + L) )
end

# ╔═╡ 000cf56d-bc91-4032-96e0-7f8d35f7b5e3
 fp_indices = combine(groupby(sfF, :datacut), srt => futurepast_index => :futurepast_index) #bar(:datacut, :futurepast_index)

# ╔═╡ 920ba85a-d5d7-420e-9512-95bbb34211fd
@df fp_indices bar(:datacut, :futurepast_index, label="", ylabel="Future skew")

# ╔═╡ fd19b849-318c-458c-8b5c-238a13ec7738


# ╔═╡ 9bd5aef9-e0be-47a4-8994-720b3eea4bb4
md"We could also clarify this curve by its sigmoidal shape. Namely fit the sigmoid and use it's anatomical features to differentiate"

# ╔═╡ 2e27c6f5-f5d2-4dbb-8842-1def4f407476
begin
	item = "Sort datasetᵢ=$(key1.datacut) by datasetⱼ=$(key2.datacut)"
md"""
# $item
"""
end

# ╔═╡ 793fe22f-0bf9-41c7-9555-1881d9d5a50b
compare=dropmissing(unstack(sfF, [:shift, :unit], :datacut, Symbol(srt))) 

# ╔═╡ 23b184e3-af82-4020-ab87-c82f69459ebd
begin
	scatter(compare[:,key1.datacut], compare[:,key2.datacut], label="")
	plot!(-2:2, -2:2, c=:white, linestyle=:dash, xlabel="$(key1.datacut)", ylabel="$(key2.datacut)", label="")
end

# ╔═╡ 4950cbeb-2d28-4429-94e7-b63cd14076e5
md"Make version where whivever axis is closer to y=x, flip that point to the x"

# ╔═╡ 220fff65-a5f0-4643-b0d8-cf67ef2c9683
begin
	L = @layout  [a 
				 [b c d e f]] 
	plot(
		Plot.timeshift.shiftedfieldplot(SF[unit], collect(-2:1:2)),
		plot(@subset(bSF, :unit .== unit).value,  label="unit $unit"),
		layout=L
	) 
end

# ╔═╡ 305dda42-217e-4bbb-b683-45e5acc016bf
SF[unit]

# ╔═╡ Cell order:
# ╟─5005402b-4d25-41a3-916b-4c814faa9065
# ╟─44079721-eb63-4c82-af86-ee3b4f55ab42
# ╟─5c008e3c-9395-4f43-be60-ff510a54b97e
# ╟─e8c8870a-0271-11ed-2ac4-317a38722303
# ╟─e4174f16-1628-4e7c-8d30-840690422582
# ╠═e6605880-9bd0-43b3-81dd-643829d32cbc
# ╠═99b0eff5-6d22-4948-91ae-5080d838580e
# ╟─ffa1cf80-e6c5-4d5f-b3fb-555018fc5f75
# ╠═01000490-e094-4c88-ac60-7b9fdeba3ccc
# ╟─c753704b-f9bc-4824-a3ec-27b6f8ff3f8b
# ╟─502701a1-9c58-48ca-b85b-b580b6fecde7
# ╟─635d367d-accb-438d-b67d-3a9e5fdd7cb7
# ╠═63772ff9-2b0c-4d30-a084-3369e3e809e0
# ╟─4cfc656c-957c-4394-833d-93af05eff81c
# ╟─d600f8bf-555b-4972-94b3-fc0a07a241c0
# ╟─21e1d7d4-0120-41b6-81c0-f98819e0dd13
# ╠═4902ca06-9f29-48c7-a3f3-8a6178ddaa4c
# ╟─58591957-b16f-4c4a-b3e7-f6c36feca331
# ╠═c1053cca-23a9-40cc-a2db-bbbb07c8a928
# ╠═1a4cc678-6900-4e07-a6dc-78073df8b06a
# ╟─fb770592-cbed-4c47-b156-e2cdde09b529
# ╟─41fbea69-ba6c-48be-8b81-244f2e01d428
# ╠═9f47b399-bc14-4369-8432-7d2b8b5e872c
# ╠═4732f7ef-cec0-4435-a2aa-370babc06cd2
# ╠═ef5b86bb-57f9-464b-9cae-acbf8933de1d
# ╠═2266bece-09f2-41a2-87dc-adfb16775361
# ╟─404416fe-70d8-456a-90e7-52de79fa2eb5
# ╟─88fa5b09-0fab-4293-89a8-02ca6ba80aa0
# ╟─231d579e-23d6-4ac7-a96d-e37366277965
# ╟─3079567e-7fb0-4269-a01e-74be600b459f
# ╠═9fda323f-46e7-436d-a2a4-8094527b1b10
# ╟─5c962c04-b8a1-46e9-8933-b65c16d76ef1
# ╠═4f98b564-6696-47f9-89cc-0c5b92137753
# ╟─81ff9d0e-6e15-4182-8dd9-ee10457348e7
# ╟─8723bf8e-93cc-4d55-a274-7699daec4044
# ╠═9ce4497d-2b82-464b-9df9-c2e568a661ac
# ╟─87171f23-3f8f-4c33-8430-2c5e13732675
# ╟─1afbdcbf-3c8d-45ba-9101-8f1f0db13512
# ╟─7ce8f6db-6696-411e-8751-06dbcf1d25da
# ╟─11f54afb-b588-49b9-81c1-8a8393084916
# ╠═53b84249-b032-4e76-9d35-6dd7856194d9
# ╟─b5f732e2-369f-4d34-b070-de4352592bec
# ╠═b509332f-cd1d-4a6d-8b79-954b02be35b4
# ╠═51e1b371-da12-403a-ba4e-ce7790bd5d2a
# ╠═a6b29b91-d0fc-46a3-8e12-47cfbb207184
# ╠═e24a55c4-2946-47a4-a12a-3ef6654865c9
# ╠═a0083bf5-a36c-4b70-a077-8758708721da
# ╠═f4baa229-883a-41df-ba42-581783029bb1
# ╠═000cf56d-bc91-4032-96e0-7f8d35f7b5e3
# ╠═920ba85a-d5d7-420e-9512-95bbb34211fd
# ╠═fd19b849-318c-458c-8b5c-238a13ec7738
# ╟─9bd5aef9-e0be-47a4-8994-720b3eea4bb4
# ╟─2e27c6f5-f5d2-4dbb-8842-1def4f407476
# ╟─793fe22f-0bf9-41c7-9555-1881d9d5a50b
# ╟─23b184e3-af82-4020-ab87-c82f69459ebd
# ╟─4950cbeb-2d28-4429-94e7-b63cd14076e5
# ╠═220fff65-a5f0-4643-b0d8-cf67ef2c9683
# ╠═305dda42-217e-4bbb-b683-45e5acc016bf
