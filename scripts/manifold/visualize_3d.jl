using Plots
import Utils
import Plot: stereoplot, setfolder
import Plot
import Utils: plotutils
using Random

using SoftGlobalScope

colorschemes = Dict(
        :cuemem => :acton,
        :correct => :PuRd_3, # for now, let's make this something clearer thoo
        :stopWell => :Dark2_6,
        :startWell => :Dark2_6
    )
n = 20_000
alpha=0.005*100_000/n
randsamp(x) = shuffle(collect(1:size(x,1)))[1:n]

# Stereoscopic views
combos = Iterators.product( (:cuemem, :correct, :stopWell, :startWell), filter(x->x.dim==3, keys(inds)))
@softscope for (bfield, key) in combos

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
    anim = @animate for i in 1:360; 
        stereoplot(eachcol(em)...; alpha, theta=i)
    end
    gif(anim, plotsdir("manifold", string(bfield), Utils.namedtup.ntopt_string(key) * ".gif"))


end
