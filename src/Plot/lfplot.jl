module lfplot

using RecipesBase
using StatsPlots
using DataFrames
using DIutils
using Infiltrator

export cycleplot
"""
    cycle_plot(lfp::DataFrame)

Checks that your lfp dataframe cycle labels are correct (what you expect)
"""
function cycleplot(lf::AbstractDataFrame; otherfield=nothing, kws=(;))
    #otherfield = otherfields isa Vector ? otherfields : [otherfields]
    otherfields = otherfield isa AbstractVector ? otherfield : [otherfield]
    kws = kws isa Vector ? kws : fill(kws, length(otherfields))
    P = @df lf[1:2500,:] begin
        Plots.plot(:time, :raw, label="raw")
    end
    if :phase in propertynames(lf)
        P = @df lf[1:2500,:] begin
            Plots.plot!(:time, mod2pi.(:phase) .+100,label="phase")
        end
    end
    P = @df lf[1:2500,:] begin
        Plots.plot!(:time, 10*:cycle, label="cycle labels")
    end
    if otherfields !== nothing 
        for (otherfield,kw) in zip(otherfields,kws)
            if otherfield === nothing; continue; end
            @assert otherfield isa Symbol "otherfield=$otherfield should be a symbol"
            @df lf[1:2500, :] begin
                Plots.plot!(:time, 100 * nannorm_extrema(lf[1:2500,otherfield],
                (-1,1)); label=string(otherfield), kw...) 
            end
        end
    end
    P
end

@userplot PlotPhaseLock
"""
    plotphaselock(spikes::DataFrame)

Carries out phase locking plot
"""
function plotphaselock(plt::PlotPhaseLock; time=:time)
    spikes = plt.args[1]
    spikes isa DataFrame ? TypeError(:spikes, "1st arg should be dataframe") : nothing
    seriestype --> :histogram
    (spikes[!,time])
end

 # Visualize
 # begin
 #     ulfp = unstack(lfp, :time, :tetrode, :raw);
 #     # ulfp = unstack(lfp, :time, :tetrode, :cycle);
 #     sort(Int.(unique(lfp.tetrode)))
 #     Plot.setfolder("lfp")
 #     @time begin
 #         s = 6_000
 #         samp = Matrix(ulfp[s:s+1_000, Not(:time)])
 #         C=cor(Matrix(ulfp[1:100_000, Not(:time)]), dims=1)
 #         xtcks = (1:size(C,1), string.(all_tetrodes))
 #         ytcks = (1:size(C,2), string.(all_tetrodes))
 #         h=heatmap(Matrix(samp), xticks=xtcks, yticks=ytcks, xrotation=90, 
 #             yrotation=0, title="Theta on PYR tetrodes")
 #         p=plot()
 #         least_corr = all_tetrodes[sortperm(mean(C,dims=2)|>vec)][1:3]
 #         most_correlated_tetrode = all_tetrodes[argmax(mean(C,dims=2)|>vec)]
 #         for (col, tet) in zip(eachcol(samp), all_tetrodes)
 #             if tet in pyr_tetrodes && tet != most_correlated_tetrode
 #                 plot!(col, label=tet,  linewidth = (tet in least_corr)  ? 
 #                     3 : 1, linestyle= (tet in least_corr) ? :solid : :dash)
 #             elseif tet == most_correlated_tetrode
 #                 plot!(col, label="M",  linewidth=5, color=:black,
 #                     linestyle=:dash)
 #             else
 #                 plot!(col, label="--", linewidth=1, color=:black,
 #                     linestyle=:dash)
 #             end
 #         end
 #         current()
 #         hc=heatmap(C, xticks=xtcks, yticks=ytcks, xrotation=90, yrotation=0)
 #         layout = @layout [[a b]; c{0.5w}]
 #         plot(h, p, hc, layout=(1,3), size=(1200,400))
 #     end
 #     Plot.save("theta_sync")
 #     ulfp = nothing; 
 #     current()
 # end

end


