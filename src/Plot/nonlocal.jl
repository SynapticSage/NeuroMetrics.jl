module nonlocal
    using DataFrames
    using StatsBase
    using StatsPlots
    using Plots
    using Infiltrator

    """
        plot_cycle_example

    yields a plot a cycle example like in Jai+Frank 2021 Fig 1
    """
    function plot_cycle_example(gs::GroupedDataFrame, gl::GroupedDataFrame, 
            cycle::DataFrameRow; cycle_horizon=(-8, 8), 
            only_cells_in_gs=true, markersize=12, spike_emphasis_color=:red, 
            lfp_emphasis_color=:red, return_stats=false)
        
        before, after = cycle_horizon
        s, l = gs[(; cycle.cycle)], gl[(; cycle.cycle)]
        S = [gs[(; cycle)] for cycle in
             collect(UnitRange((cycle.cycle .+ cycle_horizon)...))
             if (; cycle) ∈ keys(gs)]
        S = vcat(S...)
        #@infiltrate
        S = only_cells_in_gs ? 
            subset(S, :unit => c -> vec(any(c .∈ unique(s.unit)',dims=2))) : S
        L = [gl[(; cycle)] for cycle in
             collect(UnitRange((cycle.cycle .+ cycle_horizon)...))
             if (; cycle) ∈ keys(gl)]
        L = vcat(L...)

        xlabel = "time"
        xlim = minimum(L.time), maximum(L.time)
        cyclim = minimum(l.time), maximum(l.time)

        SS = @df S scatter(:time, :unit; xlim, c=:black, markersize, 
                           markershape=:vline, xticks=[], ylabel="unit", label="")
        SS = @df s scatter!(:time, :unit; xlim, c=spike_emphasis_color,
                            markersize, markershape=:vline, xticks=[],
                            ylabel="unit", label="", yformatter=y->"$(Int(round(y)))",
                            yticks=unique(s.unit))
        SS = vspan!(cyclim; c=spike_emphasis_color, alpha=0.05,label="")

        LL = @df L plot(:time, :raw, label="theta filt";
                         c=:gray, xlim, xlabel, ylabel="theta")
        LL = @df L plot!(:time, :broadraw, label="broadband"; 
                    xlim, xlabel, ylabel="theta", c="black")
        LL = @df l plot!(:time, :raw; xlim, xlabel, ylabel="theta", 
                         c=lfp_emphasis_color, label="")
        LL = vspan!(cyclim; c=lfp_emphasis_color, alpha=0.1,label="")
        #linestyle=Symbol("dash")
        #LL = @df l plot!(:time, :broadraw; xlim,xlabel,
        #                 ylabel="theta", linestyle, c=lfp_emphasis_color,label="")
        #@infiltrate

        hline!([0]; c=:darkgray, linestyle=:dash, label="", xlim)

        lay = @layout [a{0.9w}; b{0.9w, 0.7h}]
        P = plot(SS, LL, layout=lay, size=(1000, 400))

        !(return_stats) ? P :
        ( plot=P, unit=unique(S.unit), cycle=cycle.cycle)
    end
        
end
