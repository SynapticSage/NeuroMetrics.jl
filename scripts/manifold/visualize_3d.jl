using Plots
import Utils
import Plot: stereoplot, setfolder
import Plot
import Utils: plotutils
using Random
using SoftGlobalScope
using Serialization
using Munge.manifold
using ElectronDisplay
@eval Plot exts = ["png"]

desc_vars = (;feature_engineer, distance, filt)
load_manis(Main;desc_vars...)
desc = desc_manis(;desc_vars...)
qlim = [0.02, 0.96]
@info desc_vars
        for key in keys(embedding)
            inds[key] = Utils.clean.inds_quantile_filter_dims(
                                embedding[key], qlim)
        end

colorschemes = Dict(
        :cuemem => :acton,
        :correct => :PuRd_3, # for now, let's make this something clearer thoo
        :stopWell => :Dark2_6,
        :startWell => :Dark2_6
    )
n = 20_000
alpha=0.005*100_000/n
randsamp(x) = shuffle(collect(1:size(x,1)))[1:n]

if isempty(samps)
    @error "Samps empty"
end

# Stereoscopic views
combos = Iterators.product( (:cuemem, :correct, :stopWell, :startWell), filter(x->x.dim==3, keys(inds)))
@showprogress for (bfield, key) in combos
    @info (bfield, key)
    setfolder("manifold",string(bfield))
    title="key=$key\nbfield=$bfield"
    ia = randsamp(inds[key])
    #  Get embedding inside quantile bounds
    em  = embedding[key][ia,:]
    # Get behavior
    cfull =  plotutils.plotcolor(beh[:,bfield], colorschemes[bfield])
    c = cfull[samps[key.s]][ia]
    # 
    stereoplot(eachcol(em)...; c, alpha, title)
    Plot.save(key)
    # MaKE a stereoplot
    @info "Animating" bfield key
    prog = Progress(360; desc="Animation")
    anim = @animate for i in 1:360; 
        stereoplot(eachcol(em)...; alpha, theta=i)
        next!(prog)
    end
    giffile = plotsdir("manifold","manifold_$desc", string(bfield), Utils.namedtup.ntopt_string(key) * ".gif")
    mkpath(dirname(giffile))
    gif(anim, giffile)
end
