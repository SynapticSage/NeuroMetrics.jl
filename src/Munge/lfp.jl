module lfp

    using DataFrames, Statistics, DirectionalStatistics, ImageFiltering
    using DIutils.Table, DIutils, Infiltrator, MATLAB, LazyGrids, ProgressMeter, DSP
    import DIutils: Table
    import UnicodePlots
    using StatsPlots
    function __init__()
        mat"addpath(genpath('/usr/local/chronux_2_12/'))"
    end
    
    function samprate(lfp::AbstractDataFrame)
        1/median(diff(lfp.time))
    end

    # Python butterworth : TODO: replce with DSP---dsp.jl gives different results
    using PyCall
    @pyimport numpy as np
    signal = PyCall.pyimport("scipy.signal")
    butter, ff = signal.butter, signal.filtfilt
    # Define the Butterworth filter
    function butterscipy(lowcut, highcut, fs; order=8, btype="band")
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        b, a = butter(order, [low, high], btype=btype)
        return b, a
    end
    # Apply the Butterworth filter
    function butter_filter(data, lowcut, highcut, fs; order=8, btype="band")
        b, a = butterscipy(lowcut, highcut, fs; order=order, btype=btype)
        y = ff(b, a, data)
        return y
    end


    export coherence
    function coherence(lf1::DataFrame, lf2::DataFrame; average=false, 
            returndf=false)
        lf1 = groupby(lf1,:tetrode)
        lf2 = groupby(lf2,:tetrode)

        T=F=t=f=nothing
        C,phi,S12,S1,S2,tet1,tet2 = [],[],[],[],[],[],[]
        @showprogress for (l1,l2) in Iterators.product(lf1,lf2)
            t1,t2 = l1.tetrode[1], l2.tetrode[1]
            X1,X2 = l1.broadraw, l2.broadraw
            X1,X2 = X1[:,:], X2[:,:]
            mat"[$c,$p,$s12,$s1,$s2,$t,$f] =cohgramc($X1, $X2, [0.2, 0.2], struct('tapers',[2 3],'padding',-1, 'Fs', 1500, 'fpass', [1,300]))"
            push!.([tet1,tet2,phi,S1,S2,S12,C], [t1,t2,p,s1,s2,s12,c])
            T,F = t,f
        end
        cohtup = (;T,F,tet1,tet2,phi,S1,S2,S12,C)
        if average
            cohtup = _getavgcoh(;cohtup...)
        end
        if returndf
            _getdfcoh(cohtup)
        else
           cohtup 
        end
    end
    function _getavgcoh(;S1,S2,S12,C,phi,T,kws...)
        S1,S2,S12,C,phi = mean.([S1,S2,S12,C,phi])
        kws = NamedTuple(k=>v for (k,v) in zip(keys(kws),values(kws))
                   if k ∉ [:tet1, :tet2])
        (;S1,S2,S12,C,phi,T,kws...)
    end
    function _getdfcoh(;T,F,S1,S2,S12,C,phi,average=false,
            tet1=nothing,tet2=nothing,kws...)
        T, F = ndgrid(vec(T),vec(F))
        if !average
            tet1 = vcat([fill(t,size(T)) for t in tet1]...);
            tet2 = vcat([fill(t,size(T)) for t in tet2]...);
            T = vcat((T for _ in 1:length(C))...);
            F = vcat((F for _ in 1:length(C))...);
            S1  = vcat(S1...)
            S2  = vcat(S2...)
            S12 = vcat(S12...)
            phi = vcat(phi...)
            C   = vcat(C...)
            D= DataFrame(vec.([T,F,tet1,tet2]), [:time,:freq,:tet1,:tet2])
            T = F = tet1 = tet2 = nothing
            D[!,:phi] = vec(phi)
            phi = nothing
            D[!,:S1] =  vec(S1)
            S1 = nothing
            D[!,:S2] =  vec(S2)
            S2 = nothing
            S12  = nothing
            D[!,:C] =  vec(C)
            D
        else
            DataFrame(vec.([T,F,phi,S1,S2,S12,C]),
                      [:time,:freq,:phi,:S1,:S2,:S12,:C])
        end
    end

    """
        bandpass(df::DataFrame, low, high; order=5, smoothkws...)
    execute bandpass on an lfp dataframe
    # Arguments
    - `df::DataFrame`: dataframe with raw field
    - `low::Number`: low frequency cutoff
    - `high::Number`: high frequency cutoff
    - `order::Int`: order of butterworth filter
    - `smoothkws...`: keyword arguments to pass to `Table.smooth.gauss`; if
      empty, no smoothing is performed
    """
    function bandpass(df::AbstractDataFrame, low, high;
                      fs=1/median(diff(df.time)),
                      p2p::Bool=true, 
                      order=5, newname=:filt, smoothkws...)
        println("creating filter... low=$low, high=$high, order=$order")
        low, high = Float64(low), Float64(high)
        broadraw = df[!,:broadraw]
        prevtype = eltype(broadraw) |> nonmissingtype
        print("broadraw => float64...")
        broadraw =  Float64.(broadraw)
        print("filtering...")
        filt = butter_filter(broadraw, low, high, fs; order=order, btype="band")
        df[!, newname] = filt
        print("hilbert...")
        goodinds = .!ismissing.(df[!,newname])
        hilb   = DSP.hilbert(disallowmissing(convert(Vector, 
                 df[goodinds,newname])))
        print("amp and phase...")
        amp, phase = abs.(hilb), angle.(hilb)
        if p2p
            phase = t2t_to_p2p(phase)
        end
        # Find higher variance tetrodes
        to_prev(x) = prevtype <: Int ? prevtype(round(x)) : prevtype.(x)
        broadraw    = to_prev(broadraw)
        df[!,newname]           = to_prev(df[!,newname])
        df[!,newname * "amp"]   = Float32.(amp)
        df[!,newname * "phase"] = Float32.(phase)
        if !(isempty(smoothkws))
            println("smoothing...")
            df = Table.gauss(df, smoothkws...)
        end
        df
    end
    function bandpass(gdf::GroupedDataFrame, pos...; kws...)
        @showprogress for (g,df) in enumerate(gdf)
            println("filtering group $g")
            bandpass(df, pos...; kws...)
        end
    end
    function bandpass(X::AbstractVector, low::Real, high::Real; fs=1500, 
                        order=4, btype="band", p2p::Bool=true)
        low, high = Float64(low), Float64(high)
        X = convert(Vector{Float64}, X)
        filt = butter_filter(X, low, high, fs; order=order, btype=btype)
        if all(isnan, filt)
            @error "all NaNs in filter, check order=$order"
        end
        hilb = DSP.hilbert(disallowmissing(convert(Vector, filt)))
        amp, phase = abs.(hilb), angle.(hilb)
        if p2p
            phase = t2t_to_p2p(phase)
        end
        return (;filt, phase, amp)
    end

    """
        bandstop

    execute bandstop on an lfp dataframe
    """
    function bandstop(df::AbstractDataFrame, low, high; 
    fs=1/median(diff(df.time)),
    field::Symbol=:broadraw,
    newfield::Symbol="filt"*string(field),
    order=5, rounds=1, scale::Union{Nothing,Symbol}=nothing, smoothkws...)
        # Butterworth and hilbert the averaged raw (averaging introduces higher
        # frequency changes)
        low  = low isa AbstractVector ? low : [low]
        high = high isa AbstractVector ? high : [high]
        initial = Float64.(df[!,field])
        filt = butter_filter(initial, low[1], high[1], fs; order=order, 
                 btype="bandstop")
        # while (rounds -= 1) != 0
        #     for (l, h) in zip(low, high)
        #         @info "Band stopping" low=l high=h
        #         F = DSP.digitalfilter(DSP.Bandstop(l, h, fs=fs), DSP.Elliptic(5, 2.0, 100.0))
        #         initial = DSP.filtfilt(F, initial)
        #     end
        # end
        df[!,newfield] = filt
        if field == :raw
            hilb   = DSP.hilbert(df[!,field])
            df.filtamp, df.filtphase = abs.(hilb), angle.(hilb)
        end
        # Scale to the same range as the original
        if scale !== nothing
            df[!,scale] = 
            diff(collect(extrema(df[!,scale]))) .* DIutils.nannorm_extrema(df[!,field], (-1,1))
        end
        # Find higher variance tetrodes
        if !(isempty(smoothkws))
            println("smoothing... with $smoothkws")
            df = Table.gauss(df, smoothkws...)
        end
        df
    end

    """
    smooth field of an lfp dataframe
    """
    function smooth(lf::DataFrame, field; ker=5)
        ker = Kernel.gaussian((ker,))
        lf[!, "smooth" * String(field)] = imfilter(lf[!,field], ker)
        lf
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
        annotate_cycles!(lfp::GroupedDataFrame, phase_col="phase",
            method="peak-to-peak")
    """
    function annotate_cycles!(lfp::GroupedDataFrame; phase_col="phase",
                              method="peak-to-peak")
        prog = Progress(length(lfp))
        iters = enumerate(lfp) |> collect
        for (_, lf) in iters
            annotate_cycles!(lf, phase_col=phase_col, method=method)
            next!(prog)
        end
    end
    """
        annotate_cycles!(lfp::DataFrame; phase_col="phase", method="peak-to-peak")

    annotates lfp dataframe with cycle numbers
    """
    function annotate_cycles!(lfp::AbstractDataFrame; phase_col="phase", 
                              method="peak-to-peak")
        @assert issorted(lfp.time)
        phase = lfp[!, phase_col]
        lfp.phase = phase_to_radians(lfp[:,phase_col])
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
            cycle_labels = UInt32.(accumulate(+, rising_zero_point))
        lfp[!,"phase"] = mod2pi.(lfp[!,phase_col])
        elseif method == "trough-to-trough"
            step_size = median(diff(phase))
            Δₚ = [0; diff(phase)]
            falling_zero_point = [(phase[1:end-1] .>=0) .& (phase[2:end] .<0) ; false]
            #rising_zero_point = [(phase[2:end] .>=0) .& (phase[1:end-1] .<0) ; false]
            cycle_labels = UInt32.(accumulate(+, falling_zero_point))
            lfp[!,"phase′"] = mod.(lfp[!,phase_col] .- pi, 2*pi)
        else
            throw(ArgumentError("Unrecognized method=$method"))
        end
        lfp[!,"cycle"] = cycle_labels .+ 1

        return lfp
    end
    annotate_cycles(lfp::DataFrame;kws...) = annotate_cycles!(lfp;kws...)

    export annotate_cycles!
    annotate_cycles!(data::AbstractDataFrame, 
                    cycles::AbstractDataFrame)::DataFrame = 
            Table.group.annotate_periods!(data, cycles)

    """
        handle_phase_asymmetry!(lfp::AbstractDataFrame, phase_col="phase")

    remove phase asymmetry by ranking phases, and rescaling to 0,2pi,
    creating a pure sawtooth sinusoid
    """
    function handle_phase_asymmetry!(lfp::AbstractDataFrame; 
    phase_col="phase", phase_col_out="phase")
        rank = replace(lfp[!,phase_col], missing=>-Inf);
        rank = invperm(sortperm(rank));
        rank = rank ./ length(rank);
        lfp[!, phase_col_out] = rank * 2π;
        return lfp
    end

    # Visualize phase asymmetry correction
    function visualize_phase_asymmetry_correction(lfp::AbstractDataFrame)
        s = sortperm(lfp.phase); 
        resid = lfp.phase .- lfp.phase_balance;
        UnicodePlots.lineplot( disallowmissing(lfp.phase[s][1:1000:end]),
            resid[s][1:1000:end], title="Phase asymmetry correction",
            xlabel="Phase", ylabel="Residual",)
    end

    function visualize_cycles(lfp::AbstractDataFrame)
        p=@df lfp[1:2500,:] begin
            UnicodePlots.lineplot(:time, :raw, height=50, width=100)
        end
        @df lfp[1:2500,:] begin
            UnicodePlots.lineplot!(p, :time, mod2pi.(:phase) .* 30 .- 15)
            UnicodePlots.lineplot!(p, :time, 10*:cycle)
        end
        p
    end


    """
        t2t_to_p2p(phases::Vector)::Vector

    converts trought to trough phase to peak to peak phase
    """
    function t2t_to_p2p(phases::AbstractVector; zeromin::Bool=true)::Vector
        phases = phases .- minimum(phases) .+ π;
        phases = mod.(phases, 2π);
        if zeromin
            # set range [0, 2π]
            phases
        else
            # set range [-π, π]
            phases .- π;
        end
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
    # Arguments
    - `lfp::DataFrame`: lfp dataframe
    - `pos...`: position columns
    - `kws...`: keyword arguments to pass to `Table.get_periods`
    # Returns
    - `tab::DataFrame`: table with cycle start and stop times
    """
    function get_cycle_table(lfp::AbstractDataFrame, pos...; kws...)
        @assert "cycle" in names(lfp)
        if :amp ∉ names(lfp)
            println("Missing amp of theta ... taking hilbert(lfp.raw)")
            hilb = hilbert(convert(Vector{Float32},lfp.raw))
            lfp.amp = abs.(hilb)
        end
        tab = Table.get_periods(lfp, "cycle", :amp=>mean, pos...; kws...)
        return tab
    end
    """
        get_cycle_table(lfp::GroupedDataFrame, pos...; kws...)

    obtain a table with the start and stop times of each lfp cycle
    PER GROUP
    # Arguments
    - `lfp::GroupedDataFrame`: lfp dataframe
    - `pos...`: position columns
    - `kws...`: keyword arguments to pass to `Table.get_periods`
    # Returns
    - `tab::DataFrame`: table with cycle start and stop times
    """
    function get_cycle_table(lfp::GroupedDataFrame, pos...; kws...)
        TAB = Vector{DataFrame}(undef, length(lfp))
        Threads.@threads for (i,lf) in enumerate(lfp)|>collect
            TAB[i] = get_cycle_table(lf, pos...; kws...)
        end
        if length(lfp.cols) == 1
            source_col=lfp.cols[1]
            K = keys(lfp.keymap) |> collect .|> first
        else
            source_col = :source
            K = lfp.cols
        end
        return vcat(TAB...; source=source_col=>K)
    end

    """
        get_tet

    get a tetrode subset of an lfp datarfame
    """
    getTet(L::DataFrame, T::Int) = filter(:tetrode=> t->t==T, L)

end
