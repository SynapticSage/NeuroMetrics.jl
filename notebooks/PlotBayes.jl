### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ e7ca66c3-04fb-437f-9617-dbd09dcd53f9
using DrWatson; quickactivate(expanduser("~/Projects/goal-code"))
	

# ╔═╡ 9e4adc4c-7db6-400e-a0d9-ece3bc2c6ade
begin
	using JSServe
	using Mmap
	using Makie
	using WGLMakie
	Page(listen_url="0.0.0.0")
end

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

# ╔═╡ b595078d-1f22-4d3a-b4b3-499c3e9716de
begin
	
	#Fig = WGLMakie.Figure()
	#ax  = WGLMakie.Axis(Fig[1,1], xlabel="x", ylabel="y")
	#t = Observable(Int(1))
	#s = Observable(Int(20))
	#r = @lift R[$s,:,:,$t]
	#T = size(R, 4)
	#Makie.heatmap!(ax, r, colormap=(:romaO,0.8), interpolate=false)
	#Fig
	
end

# ╔═╡ de344a8b-c602-474f-a5ce-e98a1759b17a
Fig

# ╔═╡ 6f7d17f8-0e00-47cb-ad93-00a61da77861
begin
	scatter(1:4)
end

# ╔═╡ Cell order:
# ╠═e7ca66c3-04fb-437f-9617-dbd09dcd53f9
# ╠═9e4adc4c-7db6-400e-a0d9-ece3bc2c6ade
# ╠═3f2de4ff-c4a8-4809-a63f-f24fe896233b
# ╠═f24c2489-45c5-46c8-b750-a2390677c6f2
# ╠═b595078d-1f22-4d3a-b4b3-499c3e9716de
# ╠═de344a8b-c602-474f-a5ce-e98a1759b17a
# ╠═6f7d17f8-0e00-47cb-ad93-00a61da77861
# ╠═7cdf5f36-9121-4510-9339-1f09a55dd428
