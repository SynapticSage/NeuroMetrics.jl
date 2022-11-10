module nonlocal
    using DataFrames
    using StatsBase
    using StatsPlots
    using Plots
    using Infiltrator
    import Utils

    export get_cycle_example
    function get_cycle_example(gs,gl,cycle; cycle_horizon=(-8,8),
                only_cells_in_gs=true, gS=nothing)
        before, after = cycle_horizon
        C = cycle isa DataFrameRow ? cycle.cycle : cycle
        s, l = gs[(; cycle=C)], gl[(; cycle=C)]
        if gS === nothing
            S = [gs[(; cycle)] for cycle in
                 collect(UnitRange((C .+ cycle_horizon)...))
                 if (; cycle) ∈ keys(gs)]
        else
            S = [gS[(; cycle)] for cycle in
                 collect(UnitRange((C .+ cycle_horizon)...))
                 if (; cycle) ∈ keys(gs)]
        end
        S = vcat(S...)
        #@infiltrate
        S = only_cells_in_gs ? 
            subset(S, :unit => c -> vec(any(c .∈ unique(s.unit)',dims=2))) : S
        L = [gl[(; cycle)] for cycle in
             collect(UnitRange((C .+ cycle_horizon)...))
             if (; cycle) ∈ keys(gl)]
        L = vcat(L...)
        s,l,S,L
    end
    export plot_cycle_example
    """
        plot_cycle_example

    yields a plot a cycle example like in Jai+Frank 2021 Fig 1
    """
    function plot_cycle_example(gs::GroupedDataFrame, gl::GroupedDataFrame, 
            cycle::DataFrameRow; cycle_horizon=(-8, 8), 
            only_cells_in_gs=true, markersize=2, spike_emphasis_color=:red, 
            gS=nothing,
            lfp_emphasis_color=:red, return_stats=false)

        s,l,S,L = get_cycle_example(gs,gl,cycle;cycle_horizon, only_cells_in_gs, gS)
        
        xlabel = "time"
        xlim   = (minimum(L.time) - 1/1500, maximum(L.time) + 1/1500)
        cyclim = (minimum(l.time) - 1/1500, maximum(l.time) + 1/1500)
        @info cyclim

        SS = scatter(S.time, S.unit; xlim, c=:black, markersize, 
                           
                           ylabel="unit", label="")
        SS = @df s scatter!(:time, :unit; xlim, c=spike_emphasis_color,
                            markersize, markershape=:circle, xticks=[], ylim=extrema(S.unit) .* [0.9,1.1], 
                            ylabel="unit", label="", yformatter=y->"$(Int(round(y)))",
                            yticks=unique(s.unit))
        if !all(Utils.in_range(s.time, cyclim))
            @infiltrate
        end
        SS = vline!(collect(cyclim); c=spike_emphasis_color,label="")

        LL = @df L plot(:time, :raw, label="theta filt";
                         c=:gray, xlim, xlabel, ylabel="theta")
        LL = @df L plot!(:time, :broadraw, label="broadband"; 
                    xlim, xlabel, ylabel="theta", c="black")
        LL = @df l plot!(:time, :raw; xlim, xlabel, ylabel="theta", 
                         c=lfp_emphasis_color, label="")
        LL = vline!(collect(cyclim); c=lfp_emphasis_color, label="")
        #linestyle=Symbol("dash")
        #LL = @df l plot!(:time, :broadraw; xlim,xlabel,
        #                 ylabel="theta", linestyle, c=lfp_emphasis_color,label="")
        #@infiltrate

        hline!([0]; c=:darkgray, linestyle=:dash, label="", xlim)

        lay = @layout [a{0.9w}; b{0.9w, 0.7h}]
        P = plot(SS, LL, layout=lay, size=(1000, 400), link=:x)

        !(return_stats) ? P :
        ( plot=P, unit=unique(S.unit), cycle=cycle.cycle, SS=SS, LL=LL)
    end
        
end
