module lfp

    using DataFrames
    using Statistics
    using DirectionalStatistics
    using ImageFiltering
    using Table

    """
        bandpass

    execute bandpass on an lfp dataframe
    """
    function bandpass(df::DataFrame, low, high; order=5)
        # Butterworth and hilbert the averaged raw (averaging introduces higher
        # freequency changes)
        df = copy(df)
        filt = DSP.analogfilter(DSP.Bandpass(low, high, fs=1/median(diff(df.time))),
                         DSP.Butterworth(order))
        df.raw =  Float64.(df.raw)
        df.raw = DSP.filtfilt(filt.p, df.raw)
        hilb   = DSP.hilbert(df.raw)
        df.amp, df.phase = abs.(hilb), angle.(hilb)
        # Find higher variance tetrodes
        if !(isempty(kws))
            df = Table.smooth.gauss(df, kws...)
        end
        df
    end

    """
        phase_to_radians

    converts phase numbers (whichever their coordinates) into radians
    """
    function phase_to_radians(phase)
        phase = Float32.(phase)
        phase = 2*π*(phase .- minimum(phase))./diff([extrema(phase)...]) .- π
        #phase = convert.(Float16, phase)
    end

    """
        annotate_cycles

    annotates lfp dataframe with cycle numbers
    """
    function annotate_cycles(lfp::DataFrame; phase_col="phase", method="peak-to-peak")
        phase = lfp[!, phase_col]
        lfp.phase = phase_to_radians(lfp[:,"phase"])
        println("Method=$method")
        if method == "resets"
            Δₚ = [0; diff(phase)]
            reset_points = UInt32.(Δₚ .< (-1.5*π))
            cycle_labels = accumulate(+, reset_points)
        elseif method == "peak-to-peak"
            step_size = median(diff(phase))
            Δₚ = [0; diff(phase)]
            #falling_zero_point = [(phase[1:end-1] .>=0) .& (phase[2:end] .<0) ; false]
            rising_zero_point = [(phase[2:end] .>=0) .& (phase[1:end-1] .<0) ; false]
            cycle_labels = accumulate(+, rising_zero_point)
            lfp[!,"phase"] = mod2pi.(lfp[!,"phase"])
        elseif method == "trough-to-trough"
            step_size = median(diff(phase))
            Δₚ = [0; diff(phase)]
            falling_zero_point = [(phase[1:end-1] .>=0) .& (phase[2:end] .<0) ; false]
            #rising_zero_point = [(phase[2:end] .>=0) .& (phase[1:end-1] .<0) ; false]
            cycle_labels = accumulate(+, falling_zero_point)
            lfp[!,"phase′"] = mod.(lfp[!,"phase"] .- pi, 2*pi)
        else
            throw(ArgumentError("Unrecognized method=$method"))
        end
        lfp[!,"cycle"] = cycle_labels
        return lfp
    end

    """
        mean_lfp

    obtains the average LFP across a set of tetrodes
    """
    function mean_lfp(lfp; mean_fields=["phase","amp","raw"], func=Circular.median)
        lfp = groupby(lfp, :time)
        non_mean_fields = setdiff(names(lfp), mean_fields)
        lfp =combine(lfp, mean_fields.=>func.=>mean_fields, 
                     non_mean_fields.=>first.=>non_mean_fields)
        return lfp[!, Not(:tetrode)]
    end

    """
        weighted

    creates an amplitude weighted average of fields across tetrodes
    """
    function weighted_lfp(lfp::DataFrame; mean_fields=["phase","amp","raw"],
            weighting="amp")
        lfp = groupby(lfp, :time)
        non_mean_fields = setdiff(names(lfp), mean_fields)
        new = DataFrame()
        @time for lf in lfp
            item = DataFrame(lf[1,:])
            for field in mean_fields
                item[!, field] .= sum(lf[!, field] .* lf[!, weighting])/sum(lf.amp)
            end
            append!(new, item)
        end
        return new[!, Not(:tetrode)]
    end

    """
        unstack_tetrode

    unstacks tetrode dimension for an lfp dataframe
    """
    function unstack_tetrode(df::DataFrame; measure::Symbol=:phase)
        unstack(df[!, [:time, :tetrode, measure]], :tetrode, measure)
    end

    """
        get_cycle_table

    obtain a table with the start and stop times of each lfp cycle
    """
    function get_cycle_table(lfp::DataFrame, pos...; kws...)
        @assert "cycle" in names(lfp)
        tab = Table.get_periods(lfp, "cycle", :amp=>mean, pos...; kws...)
        return tab
    end

    """
        get_tet

    get a tetrode subset of an lfp datarfame
    """
    getTet(L::DataFrame, T::Int) = filter(:tetrode=> t->t==T, L)

end
