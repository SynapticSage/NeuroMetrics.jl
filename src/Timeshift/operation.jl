module operation

    export func, correct, norm_information
    using DataFrames, DataFramesMeta
    using StatsBase, StatsPlots
    using Infiltrator
    import DIutils
    import DIutils: Table
    using ProgressMeter


    ########                                                        ####### 
    ########    ,---.|        |                   |    o            #######
    ########    |   ||    ,---|    ,---.,---.,---.|--- .,---.,---.  #######
    ########    |   ||    |   |    `---.|---'|    |    ||   ||   |  #######
    ########    `---'`---'`---'    `---'`---'`---'`---'``---'`   '  #######
    ########                                                      

    """
    func

    Inputs
    ------
    main, shuffle :: DataFrame
    transform :: Pair
        Which column to apply (main,shuffle) function to and what to call the
        new column
    shuffle_is_by :: set of df columns
        which dimensions should be in the shuffle
    func_dims :: set of df columns
        dimensions that should be in the final result
        (a union of shuffle_dims and these)

    """
    function func(main, shuffle, transform::Union{Symbol,Pair}; 
            func_dims=[], shuffle_is_by=[:shuffle,:unit], metric=:info,
            stat=:median, op=.-)
            
        stat = stat      isa Symbol ? eval(timeshift, stat) : stat
        tranform = transform isa Symbol ? (transform => transform) : transform
        inputfields, outputfields = transform

        shuffle = combine(groupby(shuffle, shuffle_is_by), 
                          inputfields => stat => inputfields)

        shuffle_is_by = [b for b in shuffle_is_by if b ∈ propertynames(main)]
        by = union(shuffle_is_by, func_dims)
        main, shuffle = groupby(main, by), groupby(shuffle, by)
        @assert size(main) == size(shuffle)

        result = DataFrame()
        for (m, s) in zip(main, shuffle)
            r = m[!, Not(statopfields)]
            r[!,outputfields] = op(m[!,inputfields], s[!,inputfields])
            append!(result, r)
        end

    end

    """
    correct_signalwithshuf
    """
    function correct_signalwithshuf(Is, Ss, bonf_cells_are_sig;
            signal_frac=true, shufcorr_frac=false) 

        # Get the set of significant and insignificant cells in the shuffle
        # and signal fractions
        sig_of_interest   = subset(bonf_cells_are_sig, :sig => x -> x .== signal_frac)
        insig_of_interest = subset(bonf_cells_are_sig, :sig => x -> x .== shufcorr_frac)

        # Select out those cells from the dataframes, signal and shuffle
        C = DIutils.ismember(Is.unit, sig_of_interest.unit)
        Ic_sig   = Is[C, :]
        C = DIutils.ismember(Ss.unit, insig_of_interest.unit)
        Sc_insig = Ss[C, :]

        # Get the mean signal curve and mean correction curve
        correction_curve = combine(groupby(Sc_insig, :shift), :value=>mean)
        Ic_sig.value_cor = typeof(Ic_sig.value)(undef, size(Ic_sig,1))
        for row in eachrow(correction_curve)
            corr = Ic_sig.shift .== row.shift
            @infiltrate
            Ic_sig[corr, :value_cor] = Ic_sig[corr,:value] .- row.value_mean
        end
        sig_curve = combine(groupby(Ic_sig, :shift), 
                            :value_cor => mean,
                            :value     => mean)

        sig_curve.correct = sig_curve.value_mean - correction_curve.value_mean
        curve = @df sig_curve plot(:shift * 60, :correct, fillrange=(0,:correct), 
                           alpha=0.5, label="Corrected average info\nacross population")

        return (;Icorr=Ic_sig, curve)
    end

    """
    correct_cellwise_shuf
    """
    function correct_cellwise_shuf(Is, Ss, bonf_cells_are_sig;
            signal_frac=true, shufcorr_frac=false) 
        G = Table.group.mtg_via_commonmap(:unit, Is, Ss)
        sig_of_interest   = subset(bonf_cells_are_sig, 
                                   :sig => x -> x .== signal_frac)
        Ic_sig = DataFrame()
        for (sig, shuf) in G

            if sig.unit ∉  sig_of_interest.unit
                continue
            end
            sig.value .-= nanmean(shuf.value)
            append!(Ic_sig, sig)
        end

        sig_curve        = combine(groupby(Ic_sig,   :shift), :value => mean)
        sig_curve.correct = sig_curve.value_mean - correction_curve.value_mean
        curve = @df sig_curve plot(:shift * 60, :correct,
                                   fillrange=(0,:correct), alpha=0.5,
                                   label="Corrected average info\nacross
                                   population")
        return (;Icorr=Ic_sig, curve)
    end


    """
    Normalizes information to be between 0 and 1
Examine fields that are significant
    """
    function norm_information(Isc; value=:value, minmax=[0,1],
        removenan=true)
        @info "Normalizig by unit"
        Isc = copy(Isc)
        Δ = minmax[2] - minmax[1]
        Γ = minmax[1]
        groups = groupby(Isc, :unit)
        for group in groups
            norm(x) = Δ * ((x .- minimum(x))./(maximum(x).-minimum(x))) .+ Γ
            group[!,value] = norm(group[!,value])
        end
        Isc = combine(groups, identity)
        if removenan
            Isc = Isc[DIutils.notisnan(Isc[!,value]),:]
        else
            Isc
        end
    end

    """
    Sorts by best shift and adds best shift columns if they don't exist
    """
    function add_best_shift(I; value=:value)
        I = groupby(I, :unit)
        for unit in I
            max_ = argmax(unit[:, value])
            unit[!, :bestshift] .= unit[max_, :shift]
        end
        I = combine(I, identity)
        return I
    end
    """
             add_worst_shift(I; value=:value)
    Sorts by best shift and adds best shift columns if they don't exist
    """
    function add_worst_shift(I; value=:value)
        I = groupby(I, :unit)
        for unit in I
            max_ = argmin(unit[:, value])
            unit[!, :worstshift] .= unit[max_, :shift]
        end
        I = combine(I, identity)
        return I
    end

    function add_centroid_shift(I; value=:value)
        I = groupby(I, :unit)
        for unit in I
            shift = sum(unit[!, value].*unit.shift)./sum(unit.shift)
            unit[!, :centshift] .= shift
        end
        I = combine(I, identity)
        return I
    end

    function add_best_sig_shift(I; value=:value, sig=0.05)
        I = groupby(I, :unit)
        for unit in I
            V = unit[:,value]
            V[unit.sig .< sig] .= NaN
            if all((!).(isnan.(V)))
                replace!(V, NaN=>-Inf)
                max_ = argmax(V)
                unit[!, :bestsigshift] .= unit[max_, :shift]
            else
                unit[!, :bestsigshift] .= NaN
            end
        end
        I = combine(I, identity)
        return I
    end

    function add_stat_sig_per_unit(I; sigval=:frac, stat::Function=maximum)
        I = groupby(I, :unit)
        for unit in I
            val = stat(unit[:, sigval])
            unit[!, String(sigval) * "_" * String(Symbol(stat))] .= val
        end
        I = combine(I, identity)
        return I
    end

    function score_significance(Is, Ss; groups=[:shift,:unit])
        G = Table.group.mtg_via_commonmap(groups, Is, Ss)
        P_sig, P_non = [], []
        Ic = DataFrame()
        @showprogress for (ig, sg) in G
            cig = DataFrame(copy(ig))
            # significant
            cig[!,:frac] .= mean(cig.value .> sg.value)
            cig[!,:sig]  .= 1 .- cig.frac
            # Add to DF
            append!(Ic, cig)
        end
        return Ic
    end

    function unitfield_by_unitfield(Is, x, y)
        groups = combine(groupby(Is, :unit), first)
        s = sort(groups[:,[:unit, y, x]], y)
    end

end
