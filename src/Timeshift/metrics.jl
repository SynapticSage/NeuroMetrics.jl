module metrics

    import ..Fields: ShiftedFields, ShiftedField, DictOfShiftOfUnit
    using DataFrames

    FieldObj = Union{ShiftedField, ShiftedFields}
    TableObj = Union{DataFrame,GroupedDataFrame}

    function apply(I::DictOfShiftOfUnit, lambda::Function; kws...)
        apply(ShiftedFields(I), lambda; kws...)
    end
    function apply(I::ShiftedFields, lambda::Function; kws...)
        metrics = I.metrics
        G = groupby(metrics, :unit)
        for unit in G
            unit = lambda(unit; kws...)
        end
        I
    end
    function apply(I::ShiftedField, lambda::Function; kws...)
        metrics = I.metrics
        metrics = lambda(metrics; kws...)
        I
    end

    """
        derivative

    calculates d(rate) / d(shift) .. in other words, the
    way in which the firing field is morphing of shifts
    """
    function derivative(table::TableObj; metric)
        diffs = diff(table[:, metric])
        table[!, :bestshift] .= diffs
    end
    function derivative(field::FieldObj; metric)
        apply(field, derivative; metric)
    end

    """
        best_tau

    calculates the best shift for a shift set for a metric
    """
    function best_tau(table::TableObj; metric)
        target = Symbol("bestshift_" * String(metric))
        max_ = argmax(table[:, metric])
        table[!, target] .= table[max_, :shift]
    end
    function best_tau(field::FieldObj; metric)
        apply(field, best_tau; metric)
    end

    """
        worst_tau

    calculates the worst shift for a shift set for a metric
    """
    function worst_tau(table::TableObj; metric)
        target = Symbol("worstshift_" * String(metric))
        min_ = argmin(table[:, metric])
        table[!, target] .= table[min_, :shift]
    end
    function worst_tau(field::FieldObj; metric)
        apply(field, worst_tau; metric)
    end

    """
        best_metric

    calculates the best of a metric for a shift set
    """
    function best_metric(table::TableObj; metric)
        target = Symbol("best_" * String(metric))
        max_ = max(table[:, metric])
        table[!, target] .= max_
    end
    function best_metric(field::FieldObj; metric)
        apply(field, best_metric; metric)
    end


    """
        worst_metric

    calculates the worst of a metric for a shift set
    """
    function worst_metric(table::TableObj; metric)
        target = Symbol("worst_" * String(metric))
        min_ = min(table[:, metric])
        table[!, target] .= min_
    end
    function worst_metric(field::FieldObj; metric)
        apply(field, worst_metric; metric=metric)
    end

end

