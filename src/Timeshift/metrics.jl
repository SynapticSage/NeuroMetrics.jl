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
    function best_tau()
    end

    """
        worst_tau

    calculates the worst shift for a shift set for a metric
    """
    function worst_tau()
    end

    """
        best_metric

    calculates the best of a metric for a shift set
    """
    function best_metric()
    end


    """
        worst_metric

    calculates the worst of a metric for a shift set
    """
    function worst_metric()
    end

end

