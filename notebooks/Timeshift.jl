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

# ╔═╡ dfcfa875-d122-49f1-ab24-66c1937b3134
begin
	using Revise
	import DrWatson
	DrWatson.quickactivate(expanduser("~/Projects/goal-code"))
	push!(LOAD_PATH, DrWatson.srcdir())
end

# ╔═╡ 450738b2-3d49-4e45-9a4d-ffa1721f833a
begin
	using StatsBase
	using Term
	using REPL.TerminalMenus
	using Plots
	import Arrow
	using DataFramesMeta
	using DataFrames
	using ColorSchemes
	using StatsBase
	using Infiltrator
	"Importing packages"
end

# ╔═╡ a9b4b3d2-f318-11ec-210a-a70a7964ee72
Field, Load, Timeshift, Utils, Load, Table = begin
	using GoalFetchAnalysis
	import Timeshift
	import Field
	import Utils
	import Load
	import Table
	Field, Load, Timeshift, Utils, Load, Table
end

# ╔═╡ 435f9680-1520-468e-b97c-2ea4fb2c1ff4
using PlutoUI

# ╔═╡ 8d41c178-16ee-4881-b55c-fb80f274d7bc
PlutoUI.TableOfContents(title="Shifted field plot experiments")

# ╔═╡ 42ea762b-12ed-4eb8-ade0-3bffff593690
md"""
# Package Imports
"""

# ╔═╡ 10a1552f-ffda-4b79-8a4c-fe2864bc3ae5
md"""
Import some key GoalFetchAnalysis codes
"""

# ╔═╡ 12c97814-1d83-4f82-9f5f-891abb878e60
import Plot

# ╔═╡ bf904f0e-7387-4b73-890b-fbb5fc6137de
md"""
# Loading data
First, we're going to loadup a key dataframe of interest: **cells**
"""

# ╔═╡ cadaf555-3b90-4c5b-846b-686ce4130497
@time cells = Load.load_cells("RY16", 36)

# ╔═╡ 13b0d82f-d721-40f7-b2f6-d099d41e8897
md"""And pull in our *checkpoint* shifted field files!"""

# ╔═╡ 77039e5f-1c59-463d-8299-13dd6af9631c
md"""
Our main shifted fields info stats
"""

# ╔═╡ 7ee46c7e-eeb4-4ae2-bc4b-318322aff9fb
I = Timeshift.load_mains(dataframe=true);

# ╔═╡ 3c02c040-2fed-4bed-b985-8e78bb241455
md"""
Our shuffled shifed fields stats
"""

# ╔═╡ eb62a3a9-5ba8-4737-8214-06b3c436eac2
@time S = Timeshift.load_shuffles(dataframe=true);

# ╔═╡ edcb03f0-3e52-4fd6-adbe-48d37150ba13
md"""
Our fields at shifts
"""

# ╔═╡ 61d024da-011a-42a3-b456-19475da19e78
@time F = Timeshift.load_fields(dataframe=true);

# ╔═╡ 5bad660a-016b-4e05-83c3-38e10e3323a1
md"""
# Set the *DATACUT* and *MARGINAL*
Let's grab first let the user pick **what** to look at
"""

# ╔═╡ b553e927-d6b8-469a-90de-b1b0bf9efa11
@bind datacut PlutoUI.Radio(String.(unique(I.datacut)))

# ╔═╡ 225323c9-4ed6-42ce-987d-4d5557efaa35
@bind marginal PlutoUI.Radio(String.(unique(I.marginal)))

# ╔═╡ f77875f8-21c4-4d97-8103-cdc7d33adee3
md"""
datacut=$datacut
"""

# ╔═╡ 6399b851-1f0a-4999-ae00-e66437f5e264
md"""
marginal=$marginal
"""

# ╔═╡ d47c3d1d-394b-45bb-9d8b-1dc8ddc1b9b3
md"""
## Subset the marginals
"""

# ╔═╡ bfcf1a0a-6520-4f79-b6c2-f97d8d0f34cc
Is = transform(DataFramesMeta.@subset(I, :datacut .== Symbol(datacut), :marginal .== marginal), :shift => (x->x*60) => :shift)

# ╔═╡ ee772af0-a435-4390-9fe5-d4cd3a2aacac
Ss = transform(DataFramesMeta.@subset(S, :datacut .== Symbol(datacut), :marginal .== marginal) , :shift => (x->x*60) => :shift)

# ╔═╡ d1ae7695-1e2f-4d90-9663-37f500bfd53a
groups = [:datacut, :marginal, :unit, :shift, :area]

# ╔═╡ 07daf23c-24dc-4907-bdad-ba2c91978f8e
md"""
# Obtain signficance of fields
"""

# ╔═╡ 35ea2991-0a66-47a6-bec2-d1f44725bc8e
Is_sig = Timeshift.operation.score_significance(Is, Ss)

# ╔═╡ 131fb039-7631-429b-b327-73a40e408b59
h=histogram(Is_sig.sig, title="Pvalues of $datacut and $marginal")

# ╔═╡ 5adb69a3-ab98-401e-a344-38aada960e6d
@bind cell_sig_versus_shuffle_groups PlutoUI.Radio(["unit","unit-shift"])

# ╔═╡ 3cb16c9d-eb21-4947-8597-3991917cc7f0
cell_sig_versus_shuffle_groups_act = if cell_sig_versus_shuffle_groups == "unit"
	:unit
else
	[:unit, :shift]
end

# ╔═╡ 0713b24c-cdc4-435f-a8c1-872a071b7c50
begin
	G = Table.group.mtg_via_commonmap(cell_sig_versus_shuffle_groups_act, Is_sig, Ss)
	Gsig = filter(x->any(x[1].sig .< 0.05), G)
	Gnonsig = filter(x->all(x[1].sig .>= 0.05), G)
end;

# ╔═╡ 75a5d90c-8b21-4bd9-abd2-1c2794fe31e3
md"""
### How do significant unit bits/sp look relative to shuffle?
"""

# ╔═╡ 84a0f7ed-a1f2-4bf7-82a1-150be1acba53
@bind unit_sig PlutoUI.Slider(1:length(Gsig))

# ╔═╡ b4f43b70-1ff9-49fd-ad05-1209e4713c39
function plot_sig(unit, Gsig)
# Plot
	ig, sg = Gsig[unit][1], Gsig[unit][2]
	Cs = ecdf(sg.value)
	Ci = ecdf(ig.value)
	annotation = string(round(Float64(ig.value[1]),digits=2))
	if annotation == "NaN"
			plot()
	else
		if size(ig,1) == 1
			addon=",shift=$(ig.shift[1])"
		else
			addon=""
		end
		
	p =plot(Ci, title="unit=$(Int64(ig.unit[1])) area=$(ig.area[1]), $addon", label="Actual",
						xlabel="bits per spike", ylabel="fraction", ylim=[0,1])
	plot!(Cs, label="Shuffle")
	#plot!([ig.value)..., ig.value...], [ylims()...])
	
	#annotate!(ig.value, mean([ylims()...]), text(annotation))
	p
	end

end

# ╔═╡ 0a9df174-91a1-45fe-93c2-cb787ef6c4e9
plot_sig(unit_sig, Gsig)

# ╔═╡ ec76028d-61eb-447e-bfc6-7a2d46abd411
md"""
### What about how the non-sig cells loook?
"""

# ╔═╡ 72fdc63f-ce5d-42a7-8d51-b3bb4cbd7624
@bind unit_nonsig PlutoUI.Slider(1:length(Gnonsig))

# ╔═╡ 4539b9f8-4ba1-4576-82ca-02bc2e8151d1
plot_sig(unit_nonsig, Gnonsig)

# ╔═╡ 5ffb3ee8-141d-4ee0-8021-021c746b8bc5
md"""
# Population level shifts
"""

# ╔═╡ 4f71df4b-e203-403e-bf52-611f58f9751c
begin
	Timeshift.operation.add_best_shift(Is_sig)
	Timeshift.operation.add_best_sig_shift(Is_sig)
	Timeshift.operation.add_stat_sig_per_unit(Is_sig)
	Timeshift.operation.add_stat_sig_per_unit(Is_sig; sigval=:sig, stat=minimum)
	heatmap_unitshift = Plot.timeshift.heatmap_unitshift
	norm_info = Timeshift.operation.norm_information
end;

# ╔═╡ dc0e42c3-e684-4e16-b5a7-792f4fa71466


# ╔═╡ 3f946154-9c2a-4a31-bba3-99f3e68e0dea

plot(
	heatmap_unitshift(Is_sig;  removenonsig=true, setnanzero=true, sort_rows=:bestshift, clims=(0,15)),
    heatmap_unitshift(Is_sig; removenonsig=true, setnanzero=true, sort_rows=:bestsigshift, clims=(0,15))
)

# ╔═╡ 44f6cfdd-fd5d-4757-b497-e6b58039d95c

plot(
	heatmap_unitshift(norm_info(Is_sig);  removenonsig=true, setnanzero=true, sort_rows=:bestshift, clims=(0,1)),
    heatmap_unitshift(norm_info(Is_sig); removenonsig=true, setnanzero=true, sort_rows=:bestsigshift, clims=(0,1))
)

# ╔═╡ 9e1e2e3c-b616-4066-b58d-779b57a6a766


# ╔═╡ daf34dd4-597b-41fd-bba2-ab9795870b6f
heatmap_unitshift(normed_sig; remove_nonsig=true, sort_rows=:bestshift, dropmissingrows=true)

# ╔═╡ 33c64062-002a-4eea-8e29-e8e36784a666
begin
	jitter(n::Real, factor=0.1) = n + (.5 - rand()) * factor
	s=scatter(
		jitter.(Is_sig.bestshift, 0.1), 
		jitter.(Is_sig.bestsigshift, 0.1);
		c=[get(ColorSchemes.grayyellow, s) for s in Is_sig.unit./maximum(Is_sig.unit)],
		alpha=0.1, label="Match for all sig values",
		xlabel="Best shift", ylabel="Best significant shift");
	s
end

# ╔═╡ 07d3adea-b413-461d-a68d-0383cd5ab26b


# ╔═╡ Cell order:
# ╠═8d41c178-16ee-4881-b55c-fb80f274d7bc
# ╟─dfcfa875-d122-49f1-ab24-66c1937b3134
# ╟─42ea762b-12ed-4eb8-ade0-3bffff593690
# ╟─450738b2-3d49-4e45-9a4d-ffa1721f833a
# ╟─10a1552f-ffda-4b79-8a4c-fe2864bc3ae5
# ╟─a9b4b3d2-f318-11ec-210a-a70a7964ee72
# ╟─12c97814-1d83-4f82-9f5f-891abb878e60
# ╟─bf904f0e-7387-4b73-890b-fbb5fc6137de
# ╟─cadaf555-3b90-4c5b-846b-686ce4130497
# ╟─13b0d82f-d721-40f7-b2f6-d099d41e8897
# ╟─77039e5f-1c59-463d-8299-13dd6af9631c
# ╟─7ee46c7e-eeb4-4ae2-bc4b-318322aff9fb
# ╟─3c02c040-2fed-4bed-b985-8e78bb241455
# ╟─eb62a3a9-5ba8-4737-8214-06b3c436eac2
# ╟─edcb03f0-3e52-4fd6-adbe-48d37150ba13
# ╠═61d024da-011a-42a3-b456-19475da19e78
# ╟─5bad660a-016b-4e05-83c3-38e10e3323a1
# ╟─435f9680-1520-468e-b97c-2ea4fb2c1ff4
# ╠═b553e927-d6b8-469a-90de-b1b0bf9efa11
# ╠═225323c9-4ed6-42ce-987d-4d5557efaa35
# ╟─f77875f8-21c4-4d97-8103-cdc7d33adee3
# ╟─6399b851-1f0a-4999-ae00-e66437f5e264
# ╟─d47c3d1d-394b-45bb-9d8b-1dc8ddc1b9b3
# ╟─bfcf1a0a-6520-4f79-b6c2-f97d8d0f34cc
# ╟─ee772af0-a435-4390-9fe5-d4cd3a2aacac
# ╟─d1ae7695-1e2f-4d90-9663-37f500bfd53a
# ╟─07daf23c-24dc-4907-bdad-ba2c91978f8e
# ╟─35ea2991-0a66-47a6-bec2-d1f44725bc8e
# ╟─131fb039-7631-429b-b327-73a40e408b59
# ╟─5adb69a3-ab98-401e-a344-38aada960e6d
# ╟─3cb16c9d-eb21-4947-8597-3991917cc7f0
# ╟─0713b24c-cdc4-435f-a8c1-872a071b7c50
# ╟─75a5d90c-8b21-4bd9-abd2-1c2794fe31e3
# ╟─84a0f7ed-a1f2-4bf7-82a1-150be1acba53
# ╟─b4f43b70-1ff9-49fd-ad05-1209e4713c39
# ╟─0a9df174-91a1-45fe-93c2-cb787ef6c4e9
# ╟─ec76028d-61eb-447e-bfc6-7a2d46abd411
# ╟─72fdc63f-ce5d-42a7-8d51-b3bb4cbd7624
# ╟─4539b9f8-4ba1-4576-82ca-02bc2e8151d1
# ╟─5ffb3ee8-141d-4ee0-8021-021c746b8bc5
# ╠═4f71df4b-e203-403e-bf52-611f58f9751c
# ╠═dc0e42c3-e684-4e16-b5a7-792f4fa71466
# ╠═3f946154-9c2a-4a31-bba3-99f3e68e0dea
# ╠═44f6cfdd-fd5d-4757-b497-e6b58039d95c
# ╠═9e1e2e3c-b616-4066-b58d-779b57a6a766
# ╠═daf34dd4-597b-41fd-bba2-ab9795870b6f
# ╟─33c64062-002a-4eea-8e29-e8e36784a666
# ╟─07d3adea-b413-461d-a68d-0383cd5ab26b
