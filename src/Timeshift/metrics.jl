module metrics

    using DataFrames

    """
        derivative

    calculates d(rate) / d(shift) .. in other words, the
    way in which the firing field is morphing of shifts
    """
    function derivative()
    end

    """
        best_tau

    calculates the best shift for a shift set for a metric
    """
    function best_tau(I::ShiftedFields; metric)
        I = groupby(I, :unit)
        for unit in I
            max_ = argmax(unit[:, metric])
            unit[!, :bestshift] .= unit[max_, :shift]
        end
        I = combine(I, identity)
    end

    """
        worst_tau

    calculates the worst shift for a shift set for a metric
    """
    function worst_tau(I::ShiftedFields; metric)
        I = groupby(I, :unit)
        for unit in I
            max_ = argmin(unit[:, metric])
            unit[!, :worstshift] .= unit[max_, :shift]
        end
        I = combine(I, identity)
    end

    """
        best_metric

    calculates the best of a metric for a shift set
    """
    function best_metric(I::ShiftedFields; metric)
    end


    """
        worst_metric

    calculates the worst of a metric for a shift set
    """
    function worst_metric()
    end

end

