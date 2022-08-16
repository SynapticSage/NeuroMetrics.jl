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

# ╔═╡ e7ca66c3-04fb-437f-9617-dbd09dcd53f9
using DrWatson; quickactivate(expanduser("~/Projects/goal-code"))
	

# ╔═╡ 9e4adc4c-7db6-400e-a0d9-ece3bc2c6ade
begin
	using Mmap
end

# ╔═╡ 6765b950-aebd-4991-875c-881e5a8c2932
begin
	using Makie
	using WGLMakie
	using JSServe
	JSServe.configure_server!(
	listen_url="0.0.0.0", # or "127.0.0.1"
	)
	Page()
	nothing
end

# ╔═╡ d253fa99-0243-4458-aa07-81385f2831e1
using PlutoUI

# ╔═╡ 7cdf5f36-9121-4510-9339-1f09a55dd428
# ╠═╡ disabled = true
#=╠═╡
begin

	using Observables
	
	App() do session::Session
		n = 10
		index_slider = Slider(1:n)
		volume = rand(n, n, n)
		slice = map(index_slider) do idx
			return volume[:, :, idx]
		end
		fig = Figure()
		ax, cplot = contour(fig[1, 1], volume)
		rectplot = linesegments!(ax, Rect(-1, -1, 12, 12), linewidth=2, color=:red)
		on(index_slider) do idx
			translate!(rectplot, 0,0,idx)
		end
		heatmap(fig[1, 2], slice)
		slider = DOM.div("z-index: ", index_slider, index_slider.value)
		return JSServe.record_states(session, DOM.div(slider, fig))
	end
end
  ╠═╡ =#

# ╔═╡ 83d75b6b-9e0c-4ac9-a30c-10012de0c64c
scatter(1:4)

# ╔═╡ 3f2de4ff-c4a8-4809-a63f-f24fe896233b
function read_from_mmap(dims=4)
	mmap_filename = datadir("memmap_deathstar", "decode.bin")
    io = open(mmap_filename, "r")
    D = ([read(io, Int) for i in 1:dims]...,)  
    res = mmap(io, Array{Float16, length(D)}, D)     
    close(io)
    res
end

# ╔═╡ f24c2489-45c5-46c8-b750-a2390677c6f2
R = read_from_mmap();

# ╔═╡ b101ae60-0c05-40b5-9a1d-2b7d87026ccf
begin
	t = Observable(Int(1))
	s = Observable(Int(20))
	r = @lift R[$s,:,:,$t]
	T = size(R, 4)
end

# ╔═╡ 772e1773-20da-4fd9-ae96-1e8d2c2331c9
videoplay=false

# ╔═╡ 92d3f864-e19a-4264-92d3-42f673352796
begin
	if videoplay
		for tt = t[]:T
			t[] = tt[]
			sleep(0.02)
		end
	end
	md"Enter this cell to play as a video"
end

# ╔═╡ bce5fc69-2a09-4811-a083-167a48efd3ad
begin
	TT=@bind time PlutoUI.Slider(1:T, show_value=true)
	SS=@bind shift PlutoUI.Slider(1:size(R,1),show_value=true)
	(;TT,SS)
end

# ╔═╡ 05389194-40e2-4c73-aba0-b5a84980d028
begin
	t[] = time[]
	s[] = shift[]
end

# ╔═╡ b595078d-1f22-4d3a-b4b3-499c3e9716de
begin
	
	Fig = WGLMakie.Figure()
	ts = @lift string($t[])
	ax  = WGLMakie.Axis(Fig[1,1], xlabel="x", ylabel="y", title="Timestamp=$(ts[])")
	Makie.heatmap!(ax, r, colormap=(:romaO,0.8), interpolate=false)
	Fig
end

# ╔═╡ Cell order:
# ╠═9e4adc4c-7db6-400e-a0d9-ece3bc2c6ade
# ╠═e7ca66c3-04fb-437f-9617-dbd09dcd53f9
# ╠═6765b950-aebd-4991-875c-881e5a8c2932
# ╠═83d75b6b-9e0c-4ac9-a30c-10012de0c64c
# ╠═3f2de4ff-c4a8-4809-a63f-f24fe896233b
# ╠═f24c2489-45c5-46c8-b750-a2390677c6f2
# ╠═b101ae60-0c05-40b5-9a1d-2b7d87026ccf
# ╠═d253fa99-0243-4458-aa07-81385f2831e1
# ╟─05389194-40e2-4c73-aba0-b5a84980d028
# ╠═772e1773-20da-4fd9-ae96-1e8d2c2331c9
# ╟─92d3f864-e19a-4264-92d3-42f673352796
# ╟─bce5fc69-2a09-4811-a083-167a48efd3ad
# ╠═b595078d-1f22-4d3a-b4b3-499c3e9716de
# ╠═7cdf5f36-9121-4510-9339-1f09a55dd428
