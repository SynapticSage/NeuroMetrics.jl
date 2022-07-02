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

# ╔═╡ bd63b1b9-cf1f-498f-8254-1bd60da2682d
begin
	G = @time adaptive.ulanovsky_find_grid(beh, props; widths=3)
	O = @time adaptive.get_occupancy(beh, props, G)
end

# ╔═╡ bef016cd-26d1-4de8-a970-182fe2b92e88
# Test field abilities
multiunit = @time adaptive.get_adaptivefield(spikes, props, G, O)

# ╔═╡ 410128bf-2332-4aa2-91cb-441e0235cc4a
heatmap(multiunit.rate)

# ╔═╡ d51ce0f2-03bf-4c88-9302-9fd4bc8621eb
width_select = @bind width Slider(1:1:8, show_value=true, default=4)

# ╔═╡ b88c0ec1-b150-49be-828f-6c32bb770c48
units = adaptive.ulanovsky(spikes, beh, [:x,:y]; splitby=[:unit], widths=width)

# ╔═╡ 7c6cfeb1-2c78-4480-852b-aa06cc818f76
@bind unit PlutoUI.Slider(unique(spikes.unit))

# ╔═╡ 38cff24f-bbc1-42fd-98ae-385323c2480e
width_select

# ╔═╡ 44abcbd4-5f71-4924-b77d-9680cc96044f
plot( 
	heatmap(units[(;unit=unit)].rate, aspect_ratio=1),
	#heatmap(units[(;unit=unit)].count, aspect_ratio=1)
)

# ╔═╡ Cell order:
# ╠═44dde9e4-f9ca-11ec-1348-d968780f671c
# ╠═d33ca749-25af-4b1e-b78c-00ebf54f1327
# ╠═31082fe7-ed61-4d37-a025-77420da3f24a
# ╠═bd63b1b9-cf1f-498f-8254-1bd60da2682d
# ╠═bef016cd-26d1-4de8-a970-182fe2b92e88
# ╠═410128bf-2332-4aa2-91cb-441e0235cc4a
# ╠═d51ce0f2-03bf-4c88-9302-9fd4bc8621eb
# ╠═b88c0ec1-b150-49be-828f-6c32bb770c48
# ╠═7c6cfeb1-2c78-4480-852b-aa06cc818f76
# ╠═38cff24f-bbc1-42fd-98ae-385323c2480e
# ╠═44abcbd4-5f71-4924-b77d-9680cc96044f
