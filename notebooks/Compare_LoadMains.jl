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

      adaptive = Field.adaptive
      metrics = Field.metrics
      WIDTHS = OrderedDict(
          "x"=>2.5f0, "y"=>2.5f0, "currentPathLength"=>1f0,
          "currentAngle"=>Float32(2pi/80)
      )
end;

# ╔═╡ 5005402b-4d25-41a3-916b-4c814faa9065
md"""
### Compare Mains
Purupose: This notebook exists to do some quick comparisons between datasets. Not intended to be an analysis workhorse.
"""

# ╔═╡ 44079721-eb63-4c82-af86-ee3b4f55ab42
PlutoUI.TableOfContents()

# ╔═╡ e4174f16-1628-4e7c-8d30-840690422582
md"""
# Load
"""

# ╔═╡ 99b0eff5-6d22-4948-91ae-5080d838580e
@time I = Timeshift.load_mains();

# ╔═╡ 01000490-e094-4c88-ac60-7b9fdeba3ccc
keys(I)

# ╔═╡ c753704b-f9bc-4824-a3ec-27b6f8ff3f8b
md"""
### Keyset 1
"""

# ╔═╡ 502701a1-9c58-48ca-b85b-b580b6fecde7
begin
	keysets = collect(filter(k->k.grid .== :adaptive .&& :widths ∈ propertynames(k), keys(I)))
	keysets = string.(reorderprops(removeprops(keysets, [:first,:last,:resolution,:step, :marginal]),[:thresh, :grid, :datacut]))
	dataset_pick1 = @bind k1 PlutoUI.Radio(keysets, default=keysets[2])
end

# ╔═╡ 635d367d-accb-438d-b67d-3a9e5fdd7cb7
key1 = bestpartialmatch(keys(I), eval(Meta.parse(k1)))

# ╔═╡ 63772ff9-2b0c-4d30-a084-3369e3e809e0
# ╠═╡ show_logs = false
SF1 = to_dataframe( I[key1], key_name=:shift)

# ╔═╡ 4cfc656c-957c-4394-833d-93af05eff81c
md"""
### Keyset 2
"""

# ╔═╡ d600f8bf-555b-4972-94b3-fc0a07a241c0
begin
	keysets2 = collect(filter(k->k.grid .== :adaptive, keys(I)))
	keysets2 = string.(reorderprops(removeprops(keysets2, [:first,:last,:resolution,:step, :marginal]),[:thresh, :grid, :datacut,]))
	dataset_pick2 = @bind k2 PlutoUI.Radio(keysets, default=keysets2[2])
end

# ╔═╡ 21e1d7d4-0120-41b6-81c0-f98819e0dd13
key2 = bestpartialmatch(keys(I), eval(Meta.parse(k2)))

# ╔═╡ 4902ca06-9f29-48c7-a3f3-8a6178ddaa4c
# ╠═╡ show_logs = false
SF2 = to_dataframe( I[key2], key_name=:shift)

# ╔═╡ c1053cca-23a9-40cc-a2db-bbbb07c8a928
begin
	k1f = replace(k1, ", "=>"\n")
	k2f = replace(k2, ", "=>"\n")
end;

# ╔═╡ 404416fe-70d8-456a-90e7-52de79fa2eb5
md"""
# Examine
## Munge
"""

# ╔═╡ b2cbd98b-e847-4249-8220-885cd00c7cf2
begin
	bSF1 = sort(Timeshift.operation.add_best_shift(transform(SF1, :shift=> (x-> x .* -1) => :shift)), :unit)
	bSF2 = sort(Timeshift.operation.add_best_shift(transform(SF2, :shift=> (x-> x .* -1) => :shift)), :unit)

	nbSF1 = Timeshift.operation.norm_information(bSF1)
	nbSF2 = Timeshift.operation.norm_information(bSF2)

end;

# ╔═╡ 0d816d5a-d9a0-4a4e-80b6-0b7f127701d4
(;dataset_pick1, dataset_pick2)

# ╔═╡ 3079567e-7fb0-4269-a01e-74be600b459f
begin
	q1 = quantile(SF1.value, 0.96)
	q2 = quantile(SF2.value, 0.96)
	xlim1=(0,q1)
	xlim2=(0,q2)
	xlim3=(0,max(q1,q2))
	e1 = ecdf(SF1.value)
	e2 = ecdf(SF2.value)
	pe=plot(e1, legend=:topright, label="", xlim=xlim1)
	plot!(e2, label="", c=:yellow, xlim=xlim3)
plot(
	pe,
	histogram(SF1.value, yscale=:log10, title=k1f, xlim=xlim1),
	histogram(SF2.value, yscale=:log10, title=k2f, c=:yellow, xlim=xlim2),
	layout = @layout [a; b c]
)
end

# ╔═╡ 231d579e-23d6-4ac7-a96d-e37366277965
(;dataset_pick1, dataset_pick2)

# ╔═╡ 5c962c04-b8a1-46e9-8933-b65c16d76ef1
md"""
## Plot
"""

# ╔═╡ 7c5d56a6-f656-4ce0-8f8f-eaefbc6ff74d
md"""
Conclusions
- Doesn't change the conclusion (grid=:adaptive, width=4) vs (grid=:adaptive, width=2.5)
"""

# ╔═╡ 31b8564a-f92e-4a89-8077-03bd9cf600c5
begin
	p1 = Plot.timeshift.heatmap_unitshift(bSF1, clim=(0,quantile(bSF1.value,0.85)))
	p2 = Plot.timeshift.heatmap_unitshift(bSF2, clim=(0,quantile(bSF2.value,0.85)))
	p3 = Plot.timeshift.heatmap_unitshift(nbSF1,title="$(k1f)\ninfo normalized")
		p4 = Plot.timeshift.heatmap_unitshift(nbSF2,title="$(k2f)\n info normalized")

	plot(p1,p2,p3,p4)
end

# ╔═╡ 2e27c6f5-f5d2-4dbb-8842-1def4f407476
md"""
# Experimental plot section
"""

# ╔═╡ 23b184e3-af82-4020-ab87-c82f69459ebd
@bind unit PlutoUI.Slider(sort(unique(bSF.unit)), show_value=true)

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
# ╠═44079721-eb63-4c82-af86-ee3b4f55ab42
# ╠═e8c8870a-0271-11ed-2ac4-317a38722303
# ╟─e4174f16-1628-4e7c-8d30-840690422582
# ╠═99b0eff5-6d22-4948-91ae-5080d838580e
# ╠═01000490-e094-4c88-ac60-7b9fdeba3ccc
# ╟─c753704b-f9bc-4824-a3ec-27b6f8ff3f8b
# ╠═502701a1-9c58-48ca-b85b-b580b6fecde7
# ╠═635d367d-accb-438d-b67d-3a9e5fdd7cb7
# ╟─63772ff9-2b0c-4d30-a084-3369e3e809e0
# ╟─4cfc656c-957c-4394-833d-93af05eff81c
# ╠═d600f8bf-555b-4972-94b3-fc0a07a241c0
# ╟─21e1d7d4-0120-41b6-81c0-f98819e0dd13
# ╟─4902ca06-9f29-48c7-a3f3-8a6178ddaa4c
# ╟─c1053cca-23a9-40cc-a2db-bbbb07c8a928
# ╟─404416fe-70d8-456a-90e7-52de79fa2eb5
# ╟─b2cbd98b-e847-4249-8220-885cd00c7cf2
# ╠═0d816d5a-d9a0-4a4e-80b6-0b7f127701d4
# ╟─3079567e-7fb0-4269-a01e-74be600b459f
# ╠═231d579e-23d6-4ac7-a96d-e37366277965
# ╟─5c962c04-b8a1-46e9-8933-b65c16d76ef1
# ╟─7c5d56a6-f656-4ce0-8f8f-eaefbc6ff74d
# ╟─31b8564a-f92e-4a89-8077-03bd9cf600c5
# ╠═2e27c6f5-f5d2-4dbb-8842-1def4f407476
# ╠═23b184e3-af82-4020-ab87-c82f69459ebd
# ╠═220fff65-a5f0-4643-b0d8-cf67ef2c9683
# ╠═305dda42-217e-4bbb-b683-45e5acc016bf
