module nonlocal

    using DataFrames
    using StatsBase
    using StatsPlots
    using Plots
    using Infiltrator
    import DIutils

    export get_cycle_example
    """
        get_cycle_example

    this obtains an example theta cycle

    #TODO Issues
    occasionally renders some spikes slightly outside the cycle of interest
    11/100 cycles
    """
    function get_cycle_example(gs,gl,cycle; cycle_horizon=(-8,8),
                only_cells_in_gs=true, gS=nothing)
        #before, after = cycle_horizon
        C = cycle isa DataFrameRow ? cycle.cycle : cycle
        s, l = gs[(; cycle=C)], gl[(; cycle=C)]
        cycles = collect(UnitRange((C .+ cycle_horizon)...))

        gS = gS === nothing ? gs : gS

        S = [gs[(; cycle)] for cycle in
             cycles if (; cycle) ∈ keys(gs)]
        S = vcat(S...)
        S = only_cells_in_gs ? 
            subset(S, :unit => c -> vec(any(c .∈ unique(s.unit)',dims=2))) : S

        L = [gl[(; cycle)] for cycle in
             cycles if (; cycle) ∈ keys(gl)]
        L = vcat(L...)

        (;s,l,S,L)
    end

    function plot_cycle_stats(stats)
        lims, rel_labels = stats.lims, stats.rel_labels
        p=plot(lims.L[1,:], lims.L[1,:]; c=:black, linestyle=:dash, size=(1200,800), label="")
        vline!(lims.L[1,:]; c=:black, linestyle=:dash, label="")
        hline!(lims.L[2,:]; c=:black, linestyle=:dash, label="")
        scatter!(eachrow(lims.L)...; c=:blue, markersize=20, label="")
        txt = text.(string.(rel_labels.L),[:white])
        annotate!(eachrow(lims.L)..., txt)
        sticks!(eachrow(lims.L)...; c=:blue, label="")
        ylims!(minimum(lims.L), maximum(lims.L))
        scatter!(eachrow(lims.S)...; c=:red, markersize=8, label="")
        txt = text.(string.(rel_labels.S),[:white], 8)
        annotate!(eachrow(lims.S)..., txt)
        p
    end

    function get_cycle_stats(gs, gl, cycle; cycle_horizon=(-9, 8),
            #only_cells_in_gs=true, 
            gS=nothing)

        C = cycle isa DataFrameRow ? cycle.cycle : cycle
        s, l = gs[(; cycle=C)], gl[(; cycle=C)]
        s = (minimum(s.time), maximum(s.time))
        l = (minimum(l.time), maximum(l.time))
        cycles = collect(UnitRange((C .+ cycle_horizon)...))
        relativecycles = collect(UnitRange(cycle_horizon...))

        gS = gS === nothing ? gs : gS

        S_cycle_set = [cycle for cycle in cycles if (; cycle) ∈ keys(gs)]
        S = [(tmp=gs[(; cycle)]; [minimum(tmp.time),maximum(tmp.time)])
             for cycle in S_cycle_set]
        sa_label = [(tmp=gs[(; cycle)]; u=unique(tmp.cycle); u[1]) for cycle in
             cycles if (; cycle) ∈ keys(gs)]
        sr_label = relativecycles[DIutils.searchsortednearest.([cycles], sa_label)]
        @assert(length(S) == length(sa_label))
        
        L_cycle_set = [cycle for cycle in cycles if (; cycle) ∈ keys(gl)]
        L = [(tmp=gl[(; cycle)]; [minimum(tmp.time),maximum(tmp.time)]) 
                for cycle in L_cycle_set]
        la_label = [(tmp=gl[(; cycle)]; unique(tmp.cycle)[1]) for cycle in
             cycles if (; cycle) ∈ keys(gl)]
        lr_label = relativecycles[DIutils.searchsortednearest.([cycles], la_label)]
        @assert(length(L) == length(la_label))

        S, L = hcat(S...), hcat(L...)

        (;lims=(;s,l,S,L),
         set_labels=(;S=S_cycle_set, L=L_cycle_set),
         rel_labels=(;S=sr_label, L=sr_label),
         abs_labels=(;S=sa_label, L=la_label)
        )
    end

    export plot_cycle_example
    """
        plot_cycle_example

    yields a plot a cycle example like in Jai+Frank 2021 Fig 1
    """
    function plot_cycle_example(gs::GroupedDataFrame, gl::GroupedDataFrame, 
            cycle::DataFrameRow; cycle_horizon=(-8, 8), 
            only_cells_in_gs=true, markersize=2, spike_emphasis_color=:red, 
            gS=nothing, xlabel = "time",
            lfp_emphasis_color=:red, return_stats=false,
            plotcycstat::Bool=true)

        s,l,S,L = get_cycle_example(gs,gl,cycle; cycle_horizon, only_cells_in_gs, gS)

        xlim   = (minimum(L.time) - 1/1500, maximum(L.time) + 1/1500)
        cyclim = (minimum(l.time) - 1/1500, maximum(l.time) + 1/1500)
        ϵ = 2*mean(diff(L.time))

        range_error = !all(DIutils.in_range(s.time, cyclim .+ (-ϵ,ϵ)))
        if range_error
            @warn("spike times not within range of their corresponding cycle")
        end
        if range_error || plotcycstat
            stats = get_cycle_stats(gs,gl,cycle;cycle_horizon,gS)
            plot_cycstat = plot_cycle_stats(stats)
        else
            stats = nothing
            plot_cycstat = nothing
        end


        @info cyclim

        SS = scatter(S.time, S.unit; xlim, c=:black, markersize, 
                           ylabel="unit", label="")
        SS = @df s scatter!(:time, :unit; xlim, c=spike_emphasis_color,
                            markersize, markershape=:circle, xticks=[], ylim=extrema(S.unit) .* [0.9,1.1], 
                            ylabel="unit", label="", yformatter=y->"$(Int(round(y)))",
                            yticks=unique(s.unit))
        
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
        (;plot=P, unit=unique(S.unit), cycle=cycle.cycle, SS=SS, LL=LL,
         plotcycstat=plot_cycstat, rangeerror=range_error, stats)
    end
        
end
