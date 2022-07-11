module chrono

    using DataFrames
    using Statistics, NaNStatistics
    import Utils
    export isminutes, ensureTimescale

    function isminutes(X::DataFrame)
        Utils.dextrema(X.time)[1] < 1440.0 # assumes less than 24 hour recording
    end

    function isminutes(X::DataFrame, data::Symbol)
        if data ∈ [:beh,:behavior]
            median(diff(X.time)) < 0.015 # assume fps <= 60
        elseif data == :spikes
            isminutes(X)
        elseif datat ∈ [:lfp, :ripple, :ripples, :theta]
            isminutes(X)
        else
            @error "Not a recognized data symbol = $data"
        end
    end

    function ensureTimescale(X::DataFrame; kws...)
        if isminutes(X; kws...)
            X = transform(X, :time => (x->x.*60) => :time, copycols=false)
        else
            X
        end
    end

end
