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

# ╔═╡ e8c8870a-0271-11ed-2ac4-317a38722303
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
      using Timeshift
      import Plot     

      adaptive = Field.adaptive
      metrics = Field.metrics
      WIDTHS = OrderedDict(
          "x"=>2.5f0, "y"=>2.5f0, "currentPathLength"=>1f0,
          "currentAngle"=>Float32(2pi/80)
      )
end

# ╔═╡ 99b0eff5-6d22-4948-91ae-5080d838580e
I = Timeshift.load_mains()

# ╔═╡ 502701a1-9c58-48ca-b85b-b580b6fecde7
@bind k PlutoUI.Radio(string.(collect(
	filter(k->k.grid .== :adaptive .&& k.first.==-2.0, keys(I))
)))

# ╔═╡ cc0b81ac-f2d4-4d90-89e9-c266537df98a
key = collect(keys(I))[findall(k .== string.(collect(keys(I))))][1]

# ╔═╡ 63772ff9-2b0c-4d30-a084-3369e3e809e0
SF = ShiftedFields( I[key] )

# ╔═╡ 075199b3-7a90-4ce6-b931-f84849ca0d94
keys(I[key])

# ╔═╡ 69f9406d-7a28-4a15-a121-c4ae733a32c6
SF.metrics

# ╔═╡ 16daa369-cce3-412e-bc59-133182400b8b
Threads.nthreads()=16

# ╔═╡ 3079567e-7fb0-4269-a01e-74be600b459f
unstack(SF.metrics, :unit, :shift, :value)

# ╔═╡ b2cbd98b-e847-4249-8220-885cd00c7cf2
bSF = Timeshift.operation.add_best_shift(
	transform(SF.metrics, :shift=> (x-> x .* -1) => :shift))

# ╔═╡ 31b8564a-f92e-4a89-8077-03bd9cf600c5
Plot.timeshift.heatmap_unitshift(bSF, clim=(0,10))

# ╔═╡ 196addfa-d81e-4100-83c3-35c33e7869b6
Plot.timeshift.heatmap_unitshift(Timeshift.operation.norm_information(bSF))

# ╔═╡ Cell order:
# ╠═e8c8870a-0271-11ed-2ac4-317a38722303
# ╠═99b0eff5-6d22-4948-91ae-5080d838580e
# ╠═502701a1-9c58-48ca-b85b-b580b6fecde7
# ╠═cc0b81ac-f2d4-4d90-89e9-c266537df98a
# ╠═63772ff9-2b0c-4d30-a084-3369e3e809e0
# ╠═075199b3-7a90-4ce6-b931-f84849ca0d94
# ╠═69f9406d-7a28-4a15-a121-c4ae733a32c6
# ╠═16daa369-cce3-412e-bc59-133182400b8b
# ╠═3079567e-7fb0-4269-a01e-74be600b459f
# ╠═b2cbd98b-e847-4249-8220-885cd00c7cf2
# ╠═31b8564a-f92e-4a89-8077-03bd9cf600c5
# ╠═196addfa-d81e-4100-83c3-35c33e7869b6
