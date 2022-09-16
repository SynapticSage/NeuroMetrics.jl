"""
fieldgrad

For dealing with and getting gradients of receptive fields
"""
module fieldgrad

    using PyCall
    using Plots
    using Infiltrator
    using ImageFiltering

    export field_grad, viz_gradient, diff_field_grad

	function field_grad(f)
	    np = pyimport("numpy")
	    ff = [ff.rate for ff in f]
	    cell = cat(ff[:]...,dims=3)
	    np.gradient(cell)
	end

	function diff_field_grad(f;dims=3)
        [diff(ff; dims) for ff in field_grad(f)]
	end

    function dodiff(f;dims=3)
        [diff(ff; dims) for ff in f]
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

        if diffmode
            nabla = dodiff(nabla);
        end
	
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

		#@info size.([xx,yy,zz])

		xlim = (minimum(xx), maximum(xx)) .+ (-3, 3)
		ylim = (minimum(yy), maximum(yy)) .+ (-3, 3)

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
            U = zeros(maximum(x), maximum(y))
            V = zeros(maximum(x), maximum(y))
            ker = Kernel.gaussian((1,1))
            #@info size(U)
            @infiltrate
            for (xx,yy,uu,vv) in zip(x,y,u,v)
                U[xx,yy] = uu
                V[xx,yy] = vv
            end
            U=imfilter(U, ker)
            V=imfilter(V, ker)
            heatmap!(1:maximum(x),1:maximum(y),abs.( sqrt.(U.^2 .+ V.^2) )')
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
