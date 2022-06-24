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
	using StatsPlots
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
@time cells = Load.load_cells("RY16", 36, "*")

# ╔═╡ 13b0d82f-d721-40f7-b2f6-d099d41e8897
md"""And pull in our *checkpoint* shifted field files!"""

# ╔═╡ 77039e5f-1c59-463d-8299-13dd6af9631c
md"""
Our main shifted fields info stats
"""

# ╔═╡ 7ee46c7e-eeb4-4ae2-bc4b-318322aff9fb
_, I = Load.register(cells, Timeshift.load_mains(dataframe=true), on="unit", transfer="meanrate");

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
begin
	F = DataFrame(Timeshift.load_fields(dataframe=true))
	F.mat = Vector{Matrix}(undef, size(F,1))
	for row in eachrow(F)
		row.mat = reshape(row.value, length(row.dim_1), length(row.dim_2))
	end
	F
end;

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
Ss = transform(DataFramesMeta.@subset(S, :datacut .== Symbol(datacut), :marginal .== marginal) , :shift => (x->x*60) => :shift);

# ╔═╡ 3c5fad0a-fd7f-4766-88fc-53c2ad7bcca4
Fs = transform(DataFramesMeta.@subset(F, :datacut .== Symbol(datacut), :marginal .== marginal) , :shift => (x->x*60) => :shift);

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

# ╔═╡ 6bcc8054-4f6d-4d82-a34a-b5541235d87a
md"""
**TODO**
Sig best shift does not make that much sense

## Does bestsigshift make sense?
"""

# ╔═╡ 33c64062-002a-4eea-8e29-e8e36784a666
begin
	jitter(n::Real, factor=0.1) = n + (.5 - rand()) * factor
	s=scatter(
		jitter.(Is_sig.bestshift, 0.1), 
		jitter.(Is_sig.bestsigshift, 0.1);
		c=[get(ColorSchemes.grayyellow, s) for s in Is_sig.unit./maximum(Is_sig.unit)],
		aspect_ratio=1,
		alpha=0.1, label="Match for all sig values",
		xlabel="Best shift", ylabel="Best significant shift");
	d = @df Is_sig density(collect(Utils.skipnan(:bestsigshift .- :bestshift)), xlabel="Best sig shift - best shift")
	plot(s,d;layout=Plots.grid(2,1),label="",)
end

# ╔═╡ c8800b8a-fd8b-43b2-a1fa-bbb76879e56d
md"""
Fraction of the time that bestshift = bestsigshift: $(mean(Is_sig.bestsigshift .== Is_sig.bestshift))
"""

# ╔═╡ 2a75cc20-9bf4-4d2e-9cd5-a5d30d8d8c76
begin
	nanvals = Is_sig.bestsigshift .== NaN
	goodvals = (!).(nanvals)
	frac_of_agreement = mean(Is_sig.bestsigshift[goodvals] .== Is_sig.bestshift[goodvals])
	md"""
	Fraction of agreement for non-nan vals = $frac_of_agreement
	"""
end

# ╔═╡ 975a6d70-ae6b-4b25-8fef-34552c39c94c
md"""
# Firing rate effects

An interesting question is whether and how many effects are created by the spike rate of the neurons. This can help identify if controls are needed.
"""

# ╔═╡ b2ed7800-80f1-4d30-a00a-c76bcd8036a9
begin
	sort!(cells,:unit)
	cells_with_info=outerjoin(cells, 
		combine(groupby(Is, :unit), 
			:value => maximum => :info_best_shift), 
		on=:unit, 
		makeunique=true)
end;

# ╔═╡ 07b01f6f-f0b8-4f23-ac64-53f7dc9ea6f2
md"""
From the following, we can see when our mean rates are very low, the spatial info jumps way up rather deterministically. Therefore, we need to counter-balance rate with binsize, so our entropy doesn't increase.
"""

# ╔═╡ 1762e086-7f10-47cd-a895-fae9a772d6d5
begin
	kws_scat_ = (;aspect_ratio=1, xlim=(0,100), ylim=(0,100), label="Rate versus best shift info", xlabel="Mean Rate", ylabel="Best shift info", alpha=0.3)
	p_rate_best_outzoom = @df cells_with_info scatter(:meanrate, :info_best_shift;kws_scat_...) 
	plot!(0:30, 0:30, c=:black, linestyle=:dash, label="unity")
	p_rate_best_zoom = @df cells_with_info scatter(:meanrate, :info_best_shift;  kws_scat_..., xlim=(0,20), ylim=(0,20))
	plot!(0:30, 0:30, c=:black, linestyle=:dash, label="unity")
	plot(p_rate_best_zoom, p_rate_best_outzoom)
end

# ╔═╡ 09819531-982a-4ded-9b0c-5187aac26e97
md"""
The shift however has an effect within a given firing rate, so it's not like there's no ability to find a good shift. Rather that it increases the variance a ton in low firing rate regimes, and possibly noise dominates shift decision.
"""

# ╔═╡ d728c9d6-4571-40bb-b1df-b34d7ae0b785
begin
	norm(x) = (x.-minimum(x))./(maximum(x).-minimum(x))
	@df Is scatter(:meanrate, :shift, log2.(:value), xlim=(0,0.25), 
	c=get.([ColorSchemes.vik], norm(:shift)),
	xlabel="Mean FR", ylabel="Shift", zlabel="Log2(Information)",
		label="log₂(info)"
	)
end

# ╔═╡ b5a4f445-1cd2-476a-aeab-a362a2de3b17
Is_sig

# ╔═╡ 32ddb00f-1eac-491f-a321-3d9f48c7f70c
begin
	@df Is_sig scatter(:meanrate, log2.(:value), xlim=(0,4), 
	c=get.([ColorSchemes.vik], norm(:shift)),
	xlabel="Mean FR", ylabel="Info", zlabel="Log2(Information)",
		label="log₂(info)",
		alpha=0.8
	)
	hline!([0],c=:black,linestyle=:dash)
end

# ╔═╡ 125653c1-55de-4410-ad30-d00433008005
begin
	scat_shift_meanrate     = @df Is_sig scatter(log10.(:meanrate), :bestshift, label="", ylabel="Best Shift", xlabel="Log₁₀(mean FR)")
	scat_shift_meanrate_log = @df Is_sig scatter(:meanrate, :bestshift, label="")
	vspan!([5, xlims()[2]], c=:orange, alpha=0.3, label="Interneuron Zone", linestyle=:dash)
	plot(scat_shift_meanrate, scat_shift_meanrate_log)
end

# ╔═╡ 07d3adea-b413-461d-a68d-0383cd5ab26b


# ╔═╡ cd3925de-dbb0-4e57-9a1e-48bf8fbb109f
md"""
# Raw fields
Are field selections for shifts any good?
"""

# ╔═╡ 2bb18fdc-9080-4f0a-9d09-c2cbe0e6404a
@bind unit_field_select PlutoUI.Slider(sort(unique(Fs.unit)), show_value=true)

# ╔═╡ 398a837c-4710-4eb6-9f78-786d7173bc49
@bind shift_select PlutoUI.Slider(sort(unique(Fs.shift)), show_value=true, default=0)

# ╔═╡ 49bb383c-b3cf-448f-9906-ab1eabae3f75
md"Neuron 🧠 $unit_field_select", md"Shift 🏃 $shift_select"

# ╔═╡ 54d550c9-0e83-4cdb-bb27-88d4c6dffe83
begin
	Fi = @subset(Fs, :unit .== unit_field_select, :shift .== shift_select)[1,:]
	cell_ratemap = @subset(cells, :unit .== unit_field_select)[1,:]
	title_ratemap = "μ(Fr) = $(round(cell_ratemap.meanrate,sigdigits=2))"
	heatmap(Fi.mat; title=title_ratemap)
end

# ╔═╡ 5df33da9-bb8d-4356-99b0-a4742a20c87e
md"""
**TODO**
1. precache heatmaps
2. clim [a. per ratemap, b. const across shifts]
3. eliminate rate map some fixed distance outside of arena/home
"""

# ╔═╡ de04cae5-7de1-40c0-bcd2-761a68fd2d6d
Fs

# ╔═╡ Cell order:
# ╟─8d41c178-16ee-4881-b55c-fb80f274d7bc
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
# ╟─61d024da-011a-42a3-b456-19475da19e78
# ╟─5bad660a-016b-4e05-83c3-38e10e3323a1
# ╟─435f9680-1520-468e-b97c-2ea4fb2c1ff4
# ╠═b553e927-d6b8-469a-90de-b1b0bf9efa11
# ╠═225323c9-4ed6-42ce-987d-4d5557efaa35
# ╟─f77875f8-21c4-4d97-8103-cdc7d33adee3
# ╟─6399b851-1f0a-4999-ae00-e66437f5e264
# ╟─d47c3d1d-394b-45bb-9d8b-1dc8ddc1b9b3
# ╟─bfcf1a0a-6520-4f79-b6c2-f97d8d0f34cc
# ╠═ee772af0-a435-4390-9fe5-d4cd3a2aacac
# ╟─3c5fad0a-fd7f-4766-88fc-53c2ad7bcca4
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
# ╟─4f71df4b-e203-403e-bf52-611f58f9751c
# ╟─3f946154-9c2a-4a31-bba3-99f3e68e0dea
# ╟─44f6cfdd-fd5d-4757-b497-e6b58039d95c
# ╟─6bcc8054-4f6d-4d82-a34a-b5541235d87a
# ╠═33c64062-002a-4eea-8e29-e8e36784a666
# ╟─c8800b8a-fd8b-43b2-a1fa-bbb76879e56d
# ╟─2a75cc20-9bf4-4d2e-9cd5-a5d30d8d8c76
# ╟─975a6d70-ae6b-4b25-8fef-34552c39c94c
# ╟─b2ed7800-80f1-4d30-a00a-c76bcd8036a9
# ╟─07b01f6f-f0b8-4f23-ac64-53f7dc9ea6f2
# ╟─1762e086-7f10-47cd-a895-fae9a772d6d5
# ╟─09819531-982a-4ded-9b0c-5187aac26e97
# ╟─d728c9d6-4571-40bb-b1df-b34d7ae0b785
# ╟─b5a4f445-1cd2-476a-aeab-a362a2de3b17
# ╟─32ddb00f-1eac-491f-a321-3d9f48c7f70c
# ╟─125653c1-55de-4410-ad30-d00433008005
# ╟─07d3adea-b413-461d-a68d-0383cd5ab26b
# ╟─cd3925de-dbb0-4e57-9a1e-48bf8fbb109f
# ╟─2bb18fdc-9080-4f0a-9d09-c2cbe0e6404a
# ╟─398a837c-4710-4eb6-9f78-786d7173bc49
# ╟─49bb383c-b3cf-448f-9906-ab1eabae3f75
# ╟─54d550c9-0e83-4cdb-bb27-88d4c6dffe83
# ╟─5df33da9-bb8d-4356-99b0-a4742a20c87e
# ╠═de04cae5-7de1-40c0-bcd2-761a68fd2d6d
