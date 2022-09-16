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

# ╔═╡ b434433e-339b-11ed-3f9e-f935e0aa89a3
# ╠═╡ show_logs = false
begin
	begin
		using DrWatson
		quickactivate(expanduser("~/Projects/goal-code"))
		using GoalFetchAnalysis
	end
end

# ╔═╡ 90726dcb-2b81-47f1-bece-38de9dfbce26
using Revise

# ╔═╡ b439ebc0-339b-11ed-2c80-8b30a59b72fa
# ╠═╡ show_logs = false
begin
	using Serialization
	using Plots
	using Timeshift.shiftmetrics
	using Field.metrics
	using Timeshift
	using DimensionalData
	using Statistics, NaNStatistics, Bootstrap
	import Plot
	using StatsBase
	using PyCall
	
	F,I,shifts = deserialize(datadir("fixed_shifts_smallerbounds_higherres_2.serial"))
	keyz = Dict(key.datacut=>key for key in keys(F))
	
	parfolder="highres2"
	theme(:bright)
	function prep(f)
	    push_dims!(f)
	    pop_metric!(f, :unit)
	    push_metric!(f, metrics.bitsperspike)
	    push_metric!(f, metrics.totalcount)
	    push_metric!(f, metrics.maxrate)
	    push_shiftmetric!(f, best_tau!; metric=:bitsperspike)
	    push_metric!(f, metrics.bitsperspikeold)
	    push_shiftmetric!(f, best_tau!; metric=:bitsperspikeold)
	    f = f[vec(all(f[:totalcount] .> 50, dims=2) .&&
	              any(f[:bitsperspike] .> 0.5, dims=2)) ,:]
	end
end

# ╔═╡ 4dfe4d5e-dd4b-4f77-bf52-9bb105fcdedc
using ProgressLogging

# ╔═╡ a9edf926-71a4-4ccd-8365-c7534c5ebcf2
using Infiltrator

# ╔═╡ 0aa75cc4-1352-4688-84d8-49e5a7fa2758
using PlutoUI

# ╔═╡ a43c9736-a29c-4dc6-8d2b-4a036bdfd413
PlutoUI.TableOfContents()

# ╔═╡ d12cf45a-0bae-45d0-b37a-29bcc5386a26
md"# Helper functions 🤝"

# ╔═╡ ddd4b260-62c4-4aca-b456-acd4f99403f4
begin
	function field_grad(f)
	    np = pyimport("numpy")
	    ff = [ff.rate for ff in f]
	    cell = cat(ff[:]...,dims=3)
	    grad = np.gradient(cell)
	end

	# as: arrow head size 0-1 (fraction of arrow length;  la: arrow alpha transparency 0-1
	function arrow0!(x, y, u, v; as=0.07, lw=1, lc=:black, la=1)
	    nuv = sqrt(u^2 + v^2)
	    v1, v2 = [u;v] / nuv,  [-v;u] / nuv
	    v4 = (3*v1 + v2)/3.1623  # sqrt(10) to get unit vector
	    v5 = v4 - 2*(v4'*v2)*v2
	    v4, v5 = as*nuv*v4, as*nuv*v5
	    plot!([x,x+u], [y,y+v], lw=lw, lc=lc, la=la)
	    plot!([x+u,x+u-v5[1]], [y+v,y+v-v5[2]], lw=lw, lc=lc, la=la)
	    plot!([x+u,x+u-v4[1]], [y+v,y+v-v4[2]], lw=lw, lc=lc, la=la)
	end
	
	function viz_gradient(nabla; diffmode=false, # d/dt(∇)
								 scale=1, # how to scale an arrow
								 ck=false, # 
								 scale_byktime=false, # scale arrow by ∂t of ∇
								 skip=1, # how many frames to stride
								 original=nothing, # side with original data?
								 modrestrict=1, # restrictions to downsamp
			as=0.07, lw=1, lc=:black, la=1,  kws...)
	
	    xx,yy,zz,ii,jj,kk = ([] for _ in 1:6)
	    for (xyz, i, j, k) in zip(CartesianIndices(nabla[1]), nabla...)
	        x,y,z = xyz.I
			if (mod(x,modrestrict), mod(y,modrestrict)) == (0,0)
		        sc = scale isa Vector ? scale[z] : scale
		        push!(xx,x)
		        push!(yy,y)
		        push!(zz,z)
		        push!(ii,i*scale)
		        push!(jj,j*scale)
		        push!(kk,k*scale)
			end
	    end

		@info size.([xx,yy,zz])

		xlim = (minimum(xx), maximum(xx)) .+ (-3, 3)
		ylim = (minimum(yy), maximum(yy)) .+ (-3, 3)

		# diff mode
		if diffmode
			for (count,k) in enumerate(unique(kk))
				@exfiltrate
				if count == 1
				end
			end
		end
	
	    g1 = @gif for zi in 1:skip:maximum(zz)
	        x,y,z,i,j, k =getindex.([xx,yy,zz,ii,jj,kk], [zz.==zi])
	        #axs[zi][:quiver](x, y, i, j)
	        
			if scale_byktime
				i = i.*k
				j = j.*k
			end
			u,v = (i, j) .* scale
			q= plot( c=:black, title="$zi";
				xlim, ylim, aspect_ratio=1, legend=:none, kws...)
			arrow0!.(x,y, u, v; as, lw, lc, la)
	        #q=quiver(x,y, quiver= 0.05 .* arrw , c=:black, title="$zi";
			#	xlim, ylim, aspect_ratio=1, kws...)
			if original !== nothing
				p=plot(original[zi])
				plot(p,q; layout=grid(2,1))
			else
				q
			end
	    end
	end
	
end

# ╔═╡ b9311113-9f9c-4244-ba68-ce33aedab6ff
mod(10,2)

# ╔═╡ ad6cb75c-ae9f-46df-82af-04ee4fc9f250
X = matrixform(F[first(keys(F))]);

# ╔═╡ 00525c93-8e18-45d8-8380-38b58f7e1ffa
typeof(X)

# ╔═╡ 28fcdd87-bd1d-47ed-9321-2f43bb723dc6
push_metric!(X, metrics.bitsperspike)

# ╔═╡ c57cb684-a311-4942-b8fc-25b937e9997e
begin
	inds=sortperm([maximum(bps) for bps in eachrow(X[:bitsperspike])])
	Xs = X[inds,:];
end

# ╔═╡ 93a47ce8-fedb-448e-b50e-011e1d6af2a0
md"# A field 🏑 "

# ╔═╡ 0a80ff77-bb8a-40af-8812-16a3e9f7c40e
unit_sel = @bind un Slider(1:size(X,1), show_value=true)

# ╔═╡ 03435c8f-43c6-4de6-8967-2f43d36b343d
f=Xs[un,:];

# ╔═╡ 39e303d6-accb-4708-9be1-2beb3f71c848
plot(Xs[unit=un,shift=At(0)])

# ╔═╡ b8e6ce65-0a64-4130-9e2a-a33cb91c05a5
grad = field_grad(f);

# ╔═╡ a49594cf-1320-45ba-91e1-9174ae2490cf
unit_sel

# ╔═╡ 6a9776a9-52cb-41f4-94d3-f77fecb39ca7
revise(Munge)

# ╔═╡ fc0e619a-0fdd-46a5-b04c-d2e5685531d2
md"# Gradient ∇ₜ"

# ╔═╡ 7ff5e7d2-bafa-4884-a550-77c35ae1b2d9
plotgrad, scale, skip, az = (;plotgrad=true, scale=2, skip=1, as=0.25)

# ╔═╡ cb82b9d7-b5c7-4102-9f8b-6fd12727ec3d
qq = plotgrad ? Munge.fieldgrad.viz_gradient(grad; scale, skip, original=Xs[unit=un], modrestrict=2, lw=1, as=az, scale_byktime=false,diffmode=true) : nothing

# ╔═╡ c7759ead-f9a5-493d-aefe-e0209ade9833
plotgradset=false

# ╔═╡ 27c278d0-0ca5-459d-a52b-6a0b8da2f9c6
if plotgradset
	for i in 1:size(Xs,1)
		gq = viz_gradient(grad; scale, skip=1, original=Xs[unit=un], modrestrict=2, lw=1, as=az, scalek=true)
	end
else
	nothing
end

# ╔═╡ 84bc63f8-6cd3-4172-9547-66f002324823
md"# Diff dX/dt"

# ╔═╡ c2e8a09d-ad5f-4e8b-9841-e2ae3f75b597
serialize(datadir("test_gradient.serial"), (;Xs,))

# ╔═╡ b2e9b1ce-4874-4196-8849-47e232b9a81e
plotgraddiff=false

# ╔═╡ 0c478b09-f91f-4f11-a693-5c5f5686d346
dq = plotgraddiff ? viz_gradient(grad; scale, skip, original=Xs[unit=un], modrestrict=2, lw=1, as=az, scale_byktime=false) : nothing

# ╔═╡ 637c9cd2-ac3b-465a-9a40-8a32176765e6


# ╔═╡ Cell order:
# ╟─b434433e-339b-11ed-3f9e-f935e0aa89a3
# ╟─b439ebc0-339b-11ed-2c80-8b30a59b72fa
# ╠═4dfe4d5e-dd4b-4f77-bf52-9bb105fcdedc
# ╠═a43c9736-a29c-4dc6-8d2b-4a036bdfd413
# ╟─d12cf45a-0bae-45d0-b37a-29bcc5386a26
# ╠═a9edf926-71a4-4ccd-8365-c7534c5ebcf2
# ╠═ddd4b260-62c4-4aca-b456-acd4f99403f4
# ╠═b9311113-9f9c-4244-ba68-ce33aedab6ff
# ╠═ad6cb75c-ae9f-46df-82af-04ee4fc9f250
# ╠═00525c93-8e18-45d8-8380-38b58f7e1ffa
# ╠═0aa75cc4-1352-4688-84d8-49e5a7fa2758
# ╠═28fcdd87-bd1d-47ed-9321-2f43bb723dc6
# ╠═c57cb684-a311-4942-b8fc-25b937e9997e
# ╟─93a47ce8-fedb-448e-b50e-011e1d6af2a0
# ╟─0a80ff77-bb8a-40af-8812-16a3e9f7c40e
# ╟─03435c8f-43c6-4de6-8967-2f43d36b343d
# ╟─39e303d6-accb-4708-9be1-2beb3f71c848
# ╠═b8e6ce65-0a64-4130-9e2a-a33cb91c05a5
# ╠═a49594cf-1320-45ba-91e1-9174ae2490cf
# ╠═90726dcb-2b81-47f1-bece-38de9dfbce26
# ╠═6a9776a9-52cb-41f4-94d3-f77fecb39ca7
# ╠═fc0e619a-0fdd-46a5-b04c-d2e5685531d2
# ╠═7ff5e7d2-bafa-4884-a550-77c35ae1b2d9
# ╠═cb82b9d7-b5c7-4102-9f8b-6fd12727ec3d
# ╠═c7759ead-f9a5-493d-aefe-e0209ade9833
# ╠═27c278d0-0ca5-459d-a52b-6a0b8da2f9c6
# ╟─84bc63f8-6cd3-4172-9547-66f002324823
# ╠═c2e8a09d-ad5f-4e8b-9841-e2ae3f75b597
# ╠═b2e9b1ce-4874-4196-8849-47e232b9a81e
# ╠═0c478b09-f91f-4f11-a693-5c5f5686d346
# ╠═637c9cd2-ac3b-465a-9a40-8a32176765e6
