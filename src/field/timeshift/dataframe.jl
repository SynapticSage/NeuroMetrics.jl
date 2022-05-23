

    function to_dataframe(shifts::AbstractDict{<:Real, <:Any}; kws...) 
        table.to_dataframe(shifts; key_name="shift", kws...)
    end
    function to_dataframe(shifts::AbstractDict{<:NamedTuple, <:Any};
            kws...)
        @debug "to_dataframe"
        table.to_dataframe(shifts; kws...)
    end

    function info_to_dataframe(shifts::AbstractDict{<:Real,<:Any};
            kws...)::DataFrame
        table.to_dataframe(shifts; key_name="shift", name="info", kws...)
    end
    function info_to_dataframe(
            shifts::AbstractDict{<:Union{NamedTuple,AbstractDict},<:Any};
            kws...)::DataFrame
        table.to_dataframe(shifts; name="info", kws...)
    end
    function info_dataframe_and_cell_dataframe(place; save_cell_table="", shift_scale=:seconds, kws...)
        df = table.to_dataframe(place, key_name="shift", name="info", kws...)
        df.shift = df.shift .* -1; # beh.time-shift... negative shift is future, so correcting this
        df = if shift_scale == :minutes
            @info "minutes"
            df = transform(df, :shift => (x-> x*60) => :shift)
        elseif shift_scale != :seconds
            @error "Must be seconds or minutes"
        else
            df
        end
        df   = sort(df, [:area,:unit,:shift])
        taus = unique(sort(df.shift))
        df_imax = combine(groupby(df, [:unit, :area]), 
                          :info=>argmax, 
                          :info=>(x->taus[argmax(x)])=>:bestTau)
        df_imax = df_imax[df_imax.bestTau.!=taus[1],:] # Excluding samples with the first tau, because that's the null condition for no variation
        return df, df_imax
    end


    function fetch_best_fields(fieldInfo::DataFrame, pos...; kws...)
        beh, data = pos
        X = Dict()
        for neuron in fieldInfo.units
            D = @subset(data, :unit.==neuron)
            x = get_fields(σ(beh,fieldInfo.bestTau), D; kws...)
            push!(X,x)
        end
        return X
    end

