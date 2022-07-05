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
begin
	  using DrWatson
	  quickactivate(expanduser("~/Projects/goal-code"))
	  using GoalFetchAnalysis
	  using Plots
	  using Revise
	  props = ["x","y"]
	  adaptive = Field.adaptive
end

# ╔═╡ d33ca749-25af-4b1e-b78c-00ebf54f1327
using PlutoUI

# ╔═╡ 31082fe7-ed61-4d37-a025-77420da3f24a
begin
	
	@time spikes, beh, ripples, cells = Load.load("RY16", 36);
	@time beh, spikes  = Load.register(beh, spikes; on="time", transfer=props)
	
end

# ╔═╡ d51ce0f2-03bf-4c88-9302-9fd4bc8621eb
width_select = @bind width Slider(1:1:8, show_value=true, default=2)

# ╔═╡ efcdc2f1-5e26-4534-953e-defae4bd8603
md"""
# Grid
Let's try making a grid object
"""

# ╔═╡ 6bdcf863-9946-4ca3-ab02-fa6aebe4b91d
G = @time adaptive.find_grid(beh, props; widths=width)

# ╔═╡ 890fe951-19bb-4c9a-a905-d798bb36c57e
O = @time adaptive.get_occupancy(beh, props, G)

# ╔═╡ bef016cd-26d1-4de8-a970-182fe2b92e88
# Test field abilities
@time multiunit = @time adaptive.get_adaptivefield(spikes, props, G, O)

# ╔═╡ 410128bf-2332-4aa2-91cb-441e0235cc4a
plot(multiunit)

# ╔═╡ b88c0ec1-b150-49be-828f-6c32bb770c48
@time units = adaptive.ulanovsky(spikes, [:x,:y], G, O; splitby=[:unit], 
                                 widths=width)

# ╔═╡ 38cff24f-bbc1-42fd-98ae-385323c2480e
width_select

# ╔═╡ 7c6cfeb1-2c78-4480-852b-aa06cc818f76
@bind unit PlutoUI.Slider(unique(spikes.unit))

# ╔═╡ 44abcbd4-5f71-4924-b77d-9680cc96044f
plot(units[(;unit=unit)], aspect_ratio=1)

# ╔═╡ 2e114813-ad97-47d8-89b1-9949cd80a268
O.camerarate

# ╔═╡ Cell order:
# ╠═44dde9e4-f9ca-11ec-1348-d968780f671c
# ╠═d33ca749-25af-4b1e-b78c-00ebf54f1327
# ╟─31082fe7-ed61-4d37-a025-77420da3f24a
# ╠═d51ce0f2-03bf-4c88-9302-9fd4bc8621eb
# ╟─efcdc2f1-5e26-4534-953e-defae4bd8603
# ╠═6bdcf863-9946-4ca3-ab02-fa6aebe4b91d
# ╠═890fe951-19bb-4c9a-a905-d798bb36c57e
# ╠═bef016cd-26d1-4de8-a970-182fe2b92e88
# ╠═410128bf-2332-4aa2-91cb-441e0235cc4a
# ╠═b88c0ec1-b150-49be-828f-6c32bb770c48
# ╠═38cff24f-bbc1-42fd-98ae-385323c2480e
# ╠═7c6cfeb1-2c78-4480-852b-aa06cc818f76
# ╠═44abcbd4-5f71-4924-b77d-9680cc96044f
# ╠═2e114813-ad97-47d8-89b1-9949cd80a268
