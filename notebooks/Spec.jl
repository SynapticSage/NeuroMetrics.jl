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

# ╔═╡ a9b4b3d2-f318-11ec-210a-a70a7964ee72
Field, Load, Timeshift, Utils, Load, Table, Plot, Munge = begin
	using GoalFetchAnalysis
	import Timeshift
	import Field
	import Utils
	import Load
	import Table
	import Plot
	import Munge
	Field, Load, Timeshift, Utils, Load, Table, Plot, Munge
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
	using Dates
	using StatsPlots
	using Logging
	using FourierAnalysis
	using DSP
	using SignalAnalysis

	if Utils.in_range(hour(now()), [0,5]) ||
	   Utils.in_range(hour(now()), [20, 24])
		Plots.theme(:dark)
		theme="dark"
	else
		Plots.theme(:bright)
		theme="bright"
	end
	
	"Importing packages, theme=$theme"
end

# ╔═╡ 435f9680-1520-468e-b97c-2ea4fb2c1ff4
using PlutoUI

# ╔═╡ 8d41c178-16ee-4881-b55c-fb80f274d7bc
PlutoUI.TableOfContents(title="Non-adaptive field shifting")

# ╔═╡ 42ea762b-12ed-4eb8-ade0-3bffff593690
md"""
# Package Imports 📦
"""

# ╔═╡ bf904f0e-7387-4b73-890b-fbb5fc6137de
md"""
# Loading data  💾
First, we're going to loadup a key dataframe of interest: **cells**
"""

# ╔═╡ cadaf555-3b90-4c5b-846b-686ce4130497
begin
	@time spikes, beh, ripples, cells = Load.load("RY16", 36);
end;

# ╔═╡ 75762e88-18e5-4fce-8954-3f02b84ffd8c
begin
	lfp = Load.load_lfp("RY16", 36)
	lfp = sort(Munge.lfp.getTet(lfp, 5), [:epoch,:time])
end

# ╔═╡ c1242599-6f28-4680-ace4-13961f41a771
FourierAnalysis.blackmanSmoother()

# ╔═╡ 04e5176b-3eac-4fff-949e-40e8a1256f76
methodswith(Spectra)

# ╔═╡ a7d10b7e-6534-4017-a042-1a414a9e96bb
md"""
# Spectra 🌊
 
Obtain the spectra
"""


# ╔═╡ 0e72f2c1-5b02-4f18-8563-d30da4d87e28
md"""
## FourierAnalysis.jl
"""

# ╔═╡ fa7280f9-0515-4076-adfd-ae9e556034d0
begin
	sr, wl = 1500, length(lfp.broadraw)
	#tapering = slepians(sr, wl)
	S = spectra(lfp.broadraw, 1500, length(lfp.broadraw); smoothing=FourierAnalysis.blackmanSmoother)
end

# ╔═╡ b9bd23f9-8167-4453-b661-6ef0d3bfe108
begin
	function my_extract(S, bp)
		S = deepcopy(S)
		y = S.y[Utils.in_range(S.flabels, bp)]
		f = filter(f-> (f > bp[1]) & (f <bp[2]) , S.flabels)
		f,y
	end
		
	plot(my_extract(S, (0,50)))
end

# ╔═╡ 653b388e-aa56-4330-bf00-991518374e55


# ╔═╡ e39a5b0e-0e68-4594-b550-9b0acaa50e53
md"""
## SignalAnalysis.jl
"""

# ╔═╡ fdcbbaa3-e74f-4031-9df1-a0ba30edec75
windowing_slide = begin
	nfslide = @bind nfft PlutoUI.Slider(1:10:6000, show_value=true, default=512)
	ovslide = @bind noverlap PlutoUI.Slider(1:10:nfft, show_value=true)
	fdowns = @bind flower PlutoUI.Slider(0:100, show_value=true, default=0)
	fups = @bind fupper PlutoUI.Slider(0:100, show_value=true, default=100)
	((;ovslide,nfslide),(;fdowns, fups))
end

# ╔═╡ e451503a-5f5d-48ef-b855-bed2dedf1964
md"""
### Spectrogram
"""

# ╔═╡ c214f1d9-f8c0-420b-bb57-06e0e9b43572
begin
	y=tfd(lfp.broadraw,Spectrogram(nfft=nfft, noverlap=500, window=hanning(nfft)))
	freq =  y.freq*1500
	power = y.power[freq .> flower .&& freq .< fupper,:]
	freq = freq[freq.>flower .&& freq.<fupper]
	z = SignalAnalysis.TFD(power, freq, y.time);
end

# ╔═╡ 60e0b0da-adef-41b8-8b08-dfbeca8082b8
begin
	tslider = @bind tslide PlutoUI.Slider(1:10:length(collect(z.time))-510)
	(;windowing_slide, tslider)
end

# ╔═╡ f6253365-de6a-48bf-a51c-b8a5ed102395
begin
	times = collect(1:500) .+ tslide
	heatmap(z.time[times], z.freq, z.power[:,times])
end

# ╔═╡ 32dd7973-e199-4ca8-b9a6-4c3530e1784d
md"""
### Coherogram

Getting coherogram

How do I define a new TFD style object with a different implementation?
"""

# ╔═╡ 666d8286-79fd-47ae-a50d-fdd474fc973c


# ╔═╡ Cell order:
# ╠═8d41c178-16ee-4881-b55c-fb80f274d7bc
# ╟─dfcfa875-d122-49f1-ab24-66c1937b3134
# ╟─42ea762b-12ed-4eb8-ade0-3bffff593690
# ╟─a9b4b3d2-f318-11ec-210a-a70a7964ee72
# ╟─450738b2-3d49-4e45-9a4d-ffa1721f833a
# ╟─bf904f0e-7387-4b73-890b-fbb5fc6137de
# ╟─cadaf555-3b90-4c5b-846b-686ce4130497
# ╠═75762e88-18e5-4fce-8954-3f02b84ffd8c
# ╠═c1242599-6f28-4680-ace4-13961f41a771
# ╟─04e5176b-3eac-4fff-949e-40e8a1256f76
# ╟─a7d10b7e-6534-4017-a042-1a414a9e96bb
# ╟─0e72f2c1-5b02-4f18-8563-d30da4d87e28
# ╠═fa7280f9-0515-4076-adfd-ae9e556034d0
# ╟─435f9680-1520-468e-b97c-2ea4fb2c1ff4
# ╟─b9bd23f9-8167-4453-b661-6ef0d3bfe108
# ╠═653b388e-aa56-4330-bf00-991518374e55
# ╟─e39a5b0e-0e68-4594-b550-9b0acaa50e53
# ╟─fdcbbaa3-e74f-4031-9df1-a0ba30edec75
# ╟─e451503a-5f5d-48ef-b855-bed2dedf1964
# ╠═c214f1d9-f8c0-420b-bb57-06e0e9b43572
# ╟─60e0b0da-adef-41b8-8b08-dfbeca8082b8
# ╠═f6253365-de6a-48bf-a51c-b8a5ed102395
# ╟─32dd7973-e199-4ca8-b9a6-4c3530e1784d
# ╠═666d8286-79fd-47ae-a50d-fdd474fc973c
