module decode
using Statistics

skipnan(x) = Iterators.filter(!isnan, x)
function movingmean(dat)
    R = [rolling(mean, dat[i,j,:],5) for i in 1:size(dat,1), j in 1:size(dat,2)]
    datnew = zeros(size(dat,1),size(dat,2), length(R[1]))
    for i in 1:size(dat,1), j in 1:size(dat,2)
        datnew[i,j,:] = vec(R[i,j])
    end
    dat=datnew
end

function quantile_threshold(dat, thresh=nothing; sample_dim=3, nan_replace_val=1)
    m = [quantile(skipnan(vec(d)), thresh) 
         for d in eachslice(dat, dims=sample_dim)]
    m = reshape(m, (1,1,length(m)))
    sz = size(dat)
    dat[dat .<= m] .= NaN
    return reshape(dat, sz)
end

module raster
end
export raster

module scatter
    using Colors, ColorSchemes

    function marker_color(time::Int, T::Vector{Float64}, padding::Float64=0.2)
        time = T[time]
        in_range = (time .> (ripples.start .- padding)) .& (time .< (ripples.stop .+ padding))
        color = colorant"black"
        if any(in_range)
            area = ripples[in_range,"area"]
            if "CA1" in area
                color = color + colorant"blue"
            else
                color = color + colorant"red"
            end
        end
        return color
    end
    function marker_color(area::Any)
        color = colorant"black"
        if "CA1" in area
            color = color + colorant"blue"
        else
            color = color + colorant"red"
        end
        return color
    end
    function glow_color(time::Int, T::Vector{Float64}, padding::Float64=0.2)
        time = T[time]
        in_range = (time .> (ripples.start .- padding)) .& (time .< (ripples.stop .+ padding))
        color = colorant"white"
        if any(in_range)
            area = ripples[in_range,"area"]
            if "CA1" in area
                color = color + colorant"blue"
            else
                color = color + colorant"red"
            end
        end
        return color
    end
    function glow_color(area::Any)
        color = colorant"white"
        if "CA1" in area
            color = color + colorant"blue"
        else
            color = color + colorant"red"
        end
        return color
    end
    function glow_width(time::Int, T::Vector{Float64}, padding::Float64=0.2)
        time = T[time]
        in_range = (time .> (ripples.start .- padding)) .& (time .< (ripples.stop .+ padding))
        width = 5
        if any(in_range)
            width = Int(round(sum(ripples[in_range,"amp"])))
        end
        return width
    end
end
export scatter

module heatmap
    using ColorSchemes, Colors
    import ColorSchemeTools
    function dynamic_colormap_per_sample(basis, samples)
    end
    function static_colormap_per_sample(basis, samples; kws...)
        cmap_bases = color_per_sample(basis, samples; kws...)
        cmap_bases = [ColorSchemeTools.make_colorscheme(  x->cmap.r,
                                                          x->cmap.g,
                                                          x->cmap.b)
                      for cmap in cmap_bases]
    end
    function color_per_sample(basis, samples; kws...)
        samples = (samples .- minimum(samples)) ./ sum(extrema(samples))
        cmap_basis = colorschemes[basis]
        cmap_bases = [get(cmap_basis, sample) for sample in samples]
    end
end
export heatmap

import GR
function volume()
end
function xtk(D)
    x = D["x_position"]
    xticks = (;xticks=collect(LinRange(0,1,length(x))), xtickslabels=x)
end

end
