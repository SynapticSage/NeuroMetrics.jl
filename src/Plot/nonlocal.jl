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
            only_cells_in_gs=true, markersize=8, spike_emphasis_color=:red, 
            lfp_emphasis_color=:red, return_stats=false)
        
        before, after = cycle_horizon
        s, l = gs[(; cycle.cycle)], gl[(; cycle.cycle)]
        S = [gs[(; cycle)] for cycle in
             collect(UnitRange((cycle.cycle .+ cycle_horizon)...))
             if (; cycle) ∈ keys(gs)]
        S = vcat(S...)
        S = only_cells_in_gs ? subset(S, :unit => c -> c .∈ unique(s.unit)) : S
        L = [gl[(; cycle)] for cycle in
             collect(UnitRange((cycle.cycle .+ cycle_horizon)...))
             if (; cycle) ∈ keys(gl)]
        L = vcat(L...)

        xlabel = "time"
        xlim = minimum(L.time), maximum(L.time)

        SS = @df S scatter(:time, :unit; xlim, c=:black, markersize, markershape=:vline, xticks=[], ylabel="unit", label="")
        LL = @df L plot(:time, :raw, label="theta filt"; xlim, xlabel, ylabel="theta")
        @df L plot!(:time, :broadraw, label="broadband"; xlim, xlabel, ylabel="theta", c="black")
        plot(SS, LL, layout=lay, size=(1000, 400))

        
        plot(SS)
        SS = @df s scatter!(:time, :unit; xlim, c=spike_emphasis_color, markersize, markershape=:vline, xticks=[], ylabel="unit", label="")
        plot(LL)
        LL = @df l plot!(:time, :raw; xlim, xlabel, ylabel="theta", c=lfp_emphasis_color, label="")
        #linestyle=Symbol("dash")
        #LL = @df l plot!(:time, :broadraw; xlim,xlabel,
        #                 ylabel="theta", linestyle, c=lfp_emphasis_color,label="")

        hline!([0], c=:darkgray, linestyle=:dash, label="")

        lay = @layout [a; b{0.7h}]
        P = plot(SS, LL, layout=lay, size=(1000, 400))
        !(return_stats) ? P :
        ( plot=P, unit=unique(S.unit), cycle=cycle.cycle)
    end
        
end
