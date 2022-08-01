module shiftmetrics

    import ..Timeshift: ShiftedFields, ShiftedField, DictOfShiftOfUnit
    using DataFrames
    using Infiltrator

    FieldObj = Union{ShiftedField, ShiftedFields}
    TableObj = Union{DataFrame,GroupedDataFrame}

    export derivative, best_tau, best_metric, worst_tau, worst_metric

# =================================================================
#                      ,---.          |             |         
#                      |---|,---.,---.|    ,   .    |--- ,---.
#                      |   ||   ||   ||    |   |    |    |   |
#                      `   '|---'|---'`---'`---|    `---'`---'
#                           |    |         `---'              
#                |         o          
#      ,-.-.,---.|--- ,---..,---.,---.
#      | | ||---'|    |    ||    `---.
#      ` ' '`---'`---'`    ``---'`---'
# =================================================================
                               

    function metricapply!(shiftedfields::DictOfShiftOfUnit, lambda::Function; kws...)
        metricapply!(ShiftedFields(shiftedfields), lambda; kws...)
    end
    function metricapply!(shiftedfields::ShiftedFields, lambda::Function; kws...)
        metrics = shiftedfields.metrics
        G = groupby(metrics, :unit)
        for unit in G
            unit = lambda(unit; kws...)
        end
        shiftedfields
    end
    function metricapply!(shiftedfield::ShiftedField, lambda::Function; kws...)
        metrics = shiftedfield.metrics
        lambda(metrics; kws...)
        shiftedfield
    end

    """
        derivative

    calculates d(rate) / d(shift) .. in other words, the
    way in which the firing field is morphing of shifts
    """
    function derivative(table::TableObj; metric)
        target = Symbol(String(:deriviatve)*"_$metric")
        diffs = [0; diff(table[:, metric])]
        table[!, target] .= diffs
    end
    function derivative(field::FieldObj; metric)
        metricapply!(field, derivative; metric)
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
        metricapply!(field, best_tau; metric)
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
        metricapply!(field, worst_tau; metric)
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
        metricapply!(field, best_metric; metric)
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
        metricapply!(field, worst_metric; metric=metric)
    end

# =================================================================
#                      ,---.          |             |         
#                      |---|,---.,---.|    ,   .    |--- ,---.
#                      |   ||   ||   ||    |   |    |    |   |
#                      `   '|---'|---'`---'`---|    `---'`---'
#                           |    |         `---'              
#    ,---.o     |        |     
#    |__. .,---.|    ,---|,---.
#    |    ||---'|    |   |`---.
#    `    ``---'`---'`---'`---'
# =================================================================

    function distance_to_centroid()
    end

end
