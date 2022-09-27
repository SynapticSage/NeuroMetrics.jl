using Plots
import Utils
import Plot

# Stereoscopic views
for (area, bfield) in Iterators.product((:cuemem, :correct, :stopWell, :startWell, :trajType),
                              areas)
    cfull =  Utils.plot.plotcolor(beh[:,bfield],:vik)
    Plot.setfolder("manifold",string(bfield))

        title="area=$area, bfield=$bfield"
        ia = inds[area]
        em  = embedding[area][ia,:]
        em2 = embedding2[area][ia,:]
        c = cfull[1:nsamp][ia]

        @gif for i in 1:360; 
            stereoplot(eachcol(em)...; alpha=0.05, theta=i)
        end

        scatter(eachcol(em2)...; c, alpha=0.05, title)
        Plot.save("scatter_area=$area")

end




