module dataframe

    import Table
    export to_dataframe, info_to_dataframe
    export fetch_best_fields
    export info_dataframe_and_cell_dataframe

    default_shiftscale = :minutes #what is output by raw.load() if behavior included in output

    # ===========================
    # Dataframe preprocessing
    # ==========================

    """
    defined set of preprocessing steps. will change if I move the sign
    inversion to the actual shifting code (someday)
    """
    function _preproc(df::DataFrame; shift_scale=nothing, keep_orig::Bool=true)::DataFrame
        if shift_scale == nothing
            throw(ArgumentError("Must provide shift_scale, for now"))
        end
        df.shift = df.shift .* -1; # beh.time-shift... negative shift is future, so correcting this
        df = if shift_scale == :minutes
            @info "minutes"
            if keep_orig
                df.shiftorig = df.shift
            end
            df = transform(df, :shift => (x-> x*60) => :shift)
        elseif shift_scale != :seconds
            df
        else
            @error "Must be seconds or minutes"
        end
        sortset = setdiff([:area,:unit,:shift], propertynames(df))
        df   = sort(df, sortset)
    end

    #  ====================
    #  Types of dataframes
    #  ====================
    function to_dataframe(measurements::AbstractDict{<:Real, <:Any}; kws...) 
        _preproc(Table.to_dataframe(shifts; key_name="shift", kws...))
    end
    function to_dataframe(measurements::AbstractDict{<:NamedTuple, <:Any};
            kws...)
        @debug "to_dataframe"
        _preproc(Table.to_dataframe(measurements; kws...))
    end

    function info_to_dataframe(measurements::AbstractDict{<:Real,<:Any};
            shift_scale=nothing,
            kws...)::DataFrame
        _preproc(
                 Table.to_dataframe(measurements; 
                                    key_name="shift", name="info", kws...);
                 shift_scale=shift_scale)
    end
    function info_to_dataframe(
            measurements::AbstractDict{<:Union{NamedTuple,AbstractDict},<:Any};
            shift_scale=nothing,
            kws...)::DataFrame
        _preproc(Table.to_dataframe(measurements; name="info", kws...);
                shift_scale=shift_scale)
    end

    function info_mean(df::DataFrame)
        combine(groupby(df, [:area, :shift]), :info=>mean)
    end

    function imax(df::DataFrame; kws...)
        taus = unique(sort(df.shift))
        if :shuffle in propertynames(df)
            groups = [:unit, :area, :shuffle]
        else
            groups = [:unit, :area]
        end
        if :shiftorig ∈ propertynames(df)
            tausorig = unique(sort(df.shiftorig))
            pos = (:info=>(x->tausorig[argmax(x)]) =>:bestTauOrig,)
        else
            pos = ()
        end
        df_imax = combine(groupby(df, groups, skipmissing=true), 
                          pos...,
                          :info=>argmax, 
                          :info=>maximum, 
                          :info=>(x->taus[argmax(x)])=>:bestTau)
        df_imax = df_imax[df_imax.bestTau.!=taus[1],:] # Excluding samples with the first tau, because that's the null condition for no variation
    end
    function imax(dict::Dict; kws...)
        df = info_to_dataframe(measurements; kws...)
        df_imax(df)
    end

    # ===========================
    # Getting fields at best tau
    # ===========================

    function fetch_best_fields(fieldInfo::DataFrame, pos...; kws...)
        beh, data = pos
        X = Dict()
        for neuron in fieldInfo.units
            D = @subset(data, :unit .== neuron)
            x = get_fields(σ(beh,fieldInfo.bestTau), D; kws...)
            push!(X,x)
        end
        return X
    end

    # ===================
    # SHorcut methods
    # ==================
    """
    shortcut method
    """
    function info_dataframe_and_cell_dataframe(measurements; save_cell_table="", shift_scale=:seconds, kws...)
        df = info_to_dataframe(measurements)
        df_imax = imax(df)
        return df, df_imax
    end

    function unstack(df::DataFrame, what=:shift, measure=:info)
        DataFrames.unstack(df, what, measure)
    end

    function get_best_fields(S::DataFrame)
    end

end
