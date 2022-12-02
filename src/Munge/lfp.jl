module lfp

    using DataFrames
    using Statistics
    using DirectionalStatistics
    using ImageFiltering
    using Table
    using Infiltrator
    using MATLAB
    using LazyGrids
    using ProgressMeter

    mat"addpath(genpath('/usr/local/chronux_2_12/'))"

    export coherence
    function coherence(lf1::DataFrame, lf2::DataFrame; average=false, 
            returndf=false)
        lf1 = groupby(lf1,:tetrode)
        lf2 = groupby(lf2,:tetrode)

        T=F=t=f=nothing
        C,phi,S12,S1,S2,tet1,tet2 = [],[],[],[],[],[],[]
        @showprogress for (l1,l2) in Iterators.product(lf1,lf2)
            t1,t2 = l1.tetrode, l2.tetrode
            X1,X2 = l1.broadraw, l2.broadraw
            X1,X2 = X1[:,:], X2[:,:]
            mat"[$c,$p,$s12,$s1,$s2,$t,$f] =cohgramc($X1, $X2, [1, 0.1], struct('tapers',[2 3],'padding',-1, 'Fs', 1500, 'fpass', [1,300]))"
            @infiltrate
            push!.([tet1,tet2,phi,S1,S2,S12,C], [tet1,tet2,p,s1,s2,s12,c])
            T,F = t,f
        end
        if average
            tet1 = tet2 = [:average]
            S1,S2,S12,C,phi = mean.([S1,S2,S12,C,phi])
        end
        if returndf
            T, F = ndgrid(t,f)
            if !average
                T   = repeat(T, outer=(length(C),1))
                F   = repeat(F, outer=(length(C),1))
                S1  = vcat(S1...)
                S2  = vcat(S2...)
                S12 = vcat(S12...)
                phi = vcat(phi...)
                C   = vcat(C...)
            end
            DataFrame([T,F,tet1,tet2,phi,S1,S2,S12,C],
                      [:time,:freq,:tet1,:tet2,:phi,:S1,:S2,:S12,:C])
        else
            (;T,F,tet1,tet2,phi,S1,S2,S12,C)
        end
    end

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
            Δₚ        = [0; diff(phase)]
            #falling_zero_point = [(phase[1:end-1] .>=0) .& (phase[2:end] .<0) ; false]
            #@infiltrate
            p = phase .- median(extrema(phase))
            rising_zero_point = [(p[2:end] .>=0) .& (p[1:end-1] .<0) ; false]
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
    annotate_cycles!(lfp::DataFrame;kws...) = annotate_cycles(lfp;kws...)

    export annotate_cycles!
    annotate_cycles!(data::DataFrame, cycles::DataFrame)::DataFrame = 
            Table.group.annotate_periods!(data, cycles)
    

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
