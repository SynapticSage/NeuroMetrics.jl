"""
    notebook_comparets

compare timeshift notebook plots
"""
module notebook_compareTS

    using Plots
    using Plot: create_blank_plot
    using Timeshift.shiftmetrics
    using DataFramesMeta, DataFrames

    function plot_shift_heatmap(tmpstatmat; met, srt, normrange, key)
		plot_heatdist = heatmap(tmpstatmat.axes[2].val, 1:length(tmpstatmat.axes[1].val), 
                                tmpstatmat, ylabel="unit", 
                                title="$(key.datacut): $met by\n $srt\n", 
                                clim=(normrange ? (-Inf,Inf) : (-25,50)), 
                                colorbar_title=(normrange ? "Neuron Min(0) to Max(1)" : "Percent above median per neuron"), 
                                colorbar_titlefontrotation= 180)
		vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")

    end

    function plot_shift_histogram(sfsmet, tmpstatmat; met, srt, normrange, key)
		plot_histdist = 
			histogram(shiftdist;  xlim, normalize=:pdf, alpha=0.5,
                      label="", xlabel="shift", ylabel="pdf")
		vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
		density!(sfsmet[:,srt]; xlim, label="", ylim=(0,0.4))
        plot_histdist
    end

	function plots_sorted_mets(sfsmet, tmpstatmat; met, srt, normrange, key, 
                               annotate_pfc=true)

        plot_shift_heatmap(tmpstatmat; met, srt, normrange, key)
		
		shifts = sort(unique(sfsmet.shift))
		#@info shifts
		xlim = (minimum(shifts), maximum(shifts))
		#@info xlim
		shiftdist = dropmissing(sfsmet, srt)[!,srt]
		#@info shiftdist

		if annotate_pfc
			sfsmet_pfc = @subset(sfsmet, :area .== "PFC", 
                                 :shift .== 0, 
                                 :unit .∈ [tmpstatmat.axes[1].val])
			sort!(sfsmet_pfc, srt)
			rows_to_mark = findall(any(sfsmet_pfc.unit' .∈ tmpstatmat.axes[1].val; dims=2)[:,1])
			@assert issorted(sfsmet_pfc[:,srt])
			@assert issorted(rows_to_mark)
			@info sfsmet_pfc[:,srt]
			@assert size(rows_to_mark) == size(sfsmet_pfc[:,srt])
			annotate!(sfsmet_pfc[:,srt], rows_to_mark, text("*", :gray))
			serialize(datadir("test.serial"),(;sfsmet_pfc, rows_to_mark, tmpstatmat, srt))
		end
		

        blank = create_blank_plot()
		lay = @layout [a;b c{0.015w}]
		plot_shift_pop_overall = plot(plot_heatdist, plot_histdist, blank;
                                      layout=lay, link=:x, size=(700,700))
	end


	function side_by_side_plot(sf1, sf2; met, srt, normrange, key1, key2)
		metmatrix1 = get_metmatrix(sf1; met, srt, normrange)
		metmatrix2 = get_metmatrix(sf2; met, srt, normrange)
		plot_shift_pop_overall1, plot_shift_pop_overall2 = 
        plot_sorted_mets(sf1, metmatrix1; met, srt, normrange=true, key=key1),
        plot_sorted_mets(sf2, metmatrix2; met, srt, normrange=true, key=key2)
		plot(plot_shift_pop_overall1, plot_shift_pop_overall2, size=(800,400))
    end

    function heatmap_ts_met_by_srt()

    end

    function dist_ts_met_by_srt()
    end

end
