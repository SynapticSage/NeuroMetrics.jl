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

# ╔═╡ ee6b3fc8-084c-11ed-34fa-1f8b4ab9fe38
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
	  using Serialization
	
      adaptive = Field.adaptive
      metrics = Field.metrics
      WIDTHS = OrderedDict(
          "x"=>2.5f0, "y"=>2.5f0, "currentPathLength"=>1f0,
          "currentAngle"=>Float32(2pi/80)
      )
end;

# ╔═╡ 548340c3-378b-47df-aa4a-14a805be3df8
# ╠═╡ disabled = true
#=╠═╡
out = deserialize(datadir("exp_pro", "RY16_36_coactivity_fields.serial"))
╠═╡ =#


# ╔═╡ 98df2272-8ef6-4f38-8f0e-15b0cf0ed873
@bind k PlutoUI.Slider(1:length(keys(out)), default=1)

# ╔═╡ 68bf580f-44a9-4a2f-9e0c-fbf27b3c20a3
K = collect(keys(out))

# ╔═╡ 5e9d11d2-a4d1-4b63-bb45-61b0cd3be55e
plot(out[K[k]])

# ╔═╡ c4e4224f-7ef7-4650-b42d-84b177840e60
K[k]

# ╔═╡ 44f7daec-477b-40de-953b-e9572a24b980
md"""
Many of the coactive events look like trajectory sweeps. Some just garbage. And some tiny but strong coherent fields.

IF continue with this, will need field stats to accept and reject many
- coherence
- information
- total events
- total trajectories
- total time (number of camera frames subtended by events)
"""

# ╔═╡ Cell order:
# ╠═ee6b3fc8-084c-11ed-34fa-1f8b4ab9fe38
# ╠═548340c3-378b-47df-aa4a-14a805be3df8
# ╠═98df2272-8ef6-4f38-8f0e-15b0cf0ed873
# ╠═68bf580f-44a9-4a2f-9e0c-fbf27b3c20a3
# ╠═5e9d11d2-a4d1-4b63-bb45-61b0cd3be55e
# ╠═c4e4224f-7ef7-4650-b42d-84b177840e60
# ╠═44f7daec-477b-40de-953b-e9572a24b980
