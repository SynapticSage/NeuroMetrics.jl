module notebook_examineTS
    using Plots
    using Timeshift.shiftmetrics

	function plot_sorted_mets(sfsmet, tmpstatmat; met, srt, normrange, key)

		plot_heatdist = heatmap(tmpstatmat.axes[2].val, 1:length(tmpstatmat.axes[1].val), tmpstatmat, ylabel="unit", title="$(key.datacut): $met by\n $srt\n", clim=(normrange ? (-Inf,Inf) : (-25,50)), colorbar_title=(normrange ? "Neuron Min(0) to Max(1)" : "Percent above median per neuron"), colorbar_titlefontrotation= 180)
		vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
		
		shifts = sort(unique(sfsmet.shift))
		#@info shifts
		xlim = (minimum(shifts), maximum(shifts))
		#@info xlim
		shiftdist = dropmissing(sfsmet, srt)[!,srt]
		#@info shiftdist

		annotate_pfc = true
		if annotate_pfc
			sfsmet_pfc = @subset(sfsmet, :area .== "PFC", :shift .== 0, :unit .∈ [tmpstatmat.axes[1].val])
			sort!(sfsmet_pfc, srt)
			rows_to_mark = findall(any(sfsmet_pfc.unit' .∈ tmpstatmat.axes[1].val; dims=2)[:,1])
			@assert issorted(sfsmet_pfc[:,srt])
			@assert issorted(rows_to_mark)
			@info sfsmet_pfc[:,srt]
			@assert size(rows_to_mark) == size(sfsmet_pfc[:,srt])
			annotate!(sfsmet_pfc[:,srt], rows_to_mark, text("*", :gray))
			serialize(datadir("test.serial"),(;sfsmet_pfc, rows_to_mark, tmpstatmat, srt))
		end
		
		plot_histdist = 
			histogram(shiftdist;  xlim, normalize=:pdf, alpha=0.5, label="", xlabel="shift", ylabel="pdf")
		vline!([0], c=:black, linestyle=:dash, linewidth=2,label="")
		density!(sfsmet[:,srt]; xlim, label="", ylim=(0,0.4))

		
		blank = plot(legend=false, grid=false, framestyle=:none, background_color_inside=:match)
	
		lay = @layout [a;b c{0.015w}]
		
		plot_shift_pop_overall = plot(plot_heatdist, plot_histdist, blank, layout=lay, link=:x, size=(700,700))
		
		plot_shift_pop_overall = plot(plot_heatdist, plot_histdist, blank, layout=lay, link=:x, size=(700,700))
		plots=plot_shift_pop_overall
	end

end

