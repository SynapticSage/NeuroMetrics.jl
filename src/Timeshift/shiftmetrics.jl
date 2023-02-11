module shiftmetrics

    using DataFrames, Infiltrator, Statistics, AxisArrays, DimensionalData

    import DIutils
    import ..Timeshift: DictOfShiftOfUnit
    import ..Timeshift.types: ShiftedFields, ShiftedField

    FieldObj = Union{ShiftedField, ShiftedFields}
    TableObj = Union{DataFrame,GroupedDataFrame, AbstractDataFrame}

    export derivative!, best_tau!, best_metric!, worst_tau!, worst_metric!,
           best_tauind!
    export push_shiftmetric!

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
                               

    """
    Adding shift metrics to dimarray objects
    """
    function metricapply!(shiftedfields::DimArray, lambda::Function; metric, kws...)
        dimnum = findfirst(name(shiftedfields.dims) .== :shift)
        shifts = collect(shiftedfields.dims[dimnum])
        shiftedfields = eachslice(shiftedfields, dims=:unit)
        for shiftedfield in shiftedfields
            mets = vcat([DataFrame(Dict(:shift=>shift,metric=>r[metric]))
                         for (shift,r) in zip(shifts,shiftedfield)]...)
            n = names(mets)
            out= lambda(mets; metric)
            newname = Symbol(setdiff(names(mets), n)[1])
            setindex!.(shiftedfield, out, [newname])
        end
    end
    push_shiftmetric! = metricapply!
    function metricapply!(shiftedfields::DictOfShiftOfUnit, lambda::Function; kws...)
        metricapply!(ShiftedFields(shiftedfields), lambda; kws...)
    end
    function metricapply!(shiftedfields::ShiftedFields, lambda::Function; kws...)
        metricapply!(shiftedfields.metrics, lambda; kws...)
        shiftedfields
    end
    function metricapply!(metricDF::AbstractDataFrame, lambda::Union{Symbol,Function}; kws...)
        lambda = lambda isa Symbol ? eval(lambda) : lambda
        G = groupby(metricDF, :unit)
        for unit in G
            unit = lambda(unit; kws...)
        end
        metricDF
    end
    function metricapply!(shiftedfield::ShiftedField, lambda::Function; kws...)
        metricDF = shiftedfield.metrics
        lambda(metricDF; kws...)
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
    function best_tau!(table::TableObj; metric)
        target = Symbol("bestshift_" * String(metric))
        max_ = argmax(table[:, metric])
        table[!, target] .= table[max_, :shift]
    end
    function best_tau!(field::FieldObj; metric)
        metricapply!(field, best_tau!; metric)
    end

    """
        best_tauind

    calculates the best shift for a shift set for a metric
    """
    function best_tauind!(table::TableObj; metric)
        target = Symbol("bestshiftind_" * String(metric))
        max_ = argmax(table[:, metric])
        table[!, target] .= max_
    end
    function best_tauind!(field::FieldObj; metric)
        metricapply!(field, best_tau!; metric)
    end

    """
        worst_tau

    calculates the worst shift for a shift set for a metric
    """
    function worst_tau!(table::TableObj; metric)
        target = Symbol("worstshift_" * String(metric))
        min_ = argmin(table[:, metric])
        table[!, target] .= table[min_, :shift]
    end
    function worst_tau!(field::FieldObj; metric)
        metricapply!(field, worst_tau!; metric)
    end

    """
        best_metric

    calculates the best of a metric for a shift set
    """
    function best_metric!(table::TableObj; metric)
        target = Symbol("best_" * String(metric))
        max_ = maximum(table[:, metric])
        table[!, target] .= max_
    end
    function best_metric!(field::FieldObj; metric)
        metricapply!(field, best_metric!; metric)
    end


    """
        worst_metric

    calculates the worst of a metric for a shift set
    """
    function worst_metric!(table::TableObj; metric)
        target = Symbol("worst_" * String(metric))
        min_ = minimum(table[:, metric])
        table[!, target] .= min_
    end
    function worst_metric!(field::FieldObj; metric)
        metricapply!(field, worst_metric!; metric=metric)
    end


    function centroid_displacement!(table::TableObj; metric=:bitsperspike)
        target = Symbol("centdisplace_" * String(metric))
        sample = Symbol("bestshift_" * String(metric))
        if sample ∉ propertynames(table)
            best_tau!(table; metric)
        end
        bestshift = table.shift .== table[!,sample][1]
        zeroshift = table.shift .== 0
        x      = table[bestshift, "centroid"] - table[zeroshift, "centroid"]
        table[!, target] .= x
        nothing
    end
    function centroid_displacement!(field::FieldObj; metric=:bitsperspike)
        metricapply!(field, centroid_displacement!; metric)
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

    getunitstat_table(sf::FieldObj) = combine(groupby(sf.metrics,
                                         intersect(propertynames(sf.metrics),
                                                   [:unit])), 
                                              x->x[x.shift.==0, :])

    function getstatmat(shiftedfields::DictOfShiftOfUnit, stat; kws...)
        getstatmat(ShiftedFields(shiftedfields), stat; kws...)
    end

    function getstatmat(shiftedfields::ShiftedFields, stat; kws...)
        getstatmat(shiftedfields, stat; kws...)
    end

    function getstatmat(sfm::DataFrame, stat;
            filtval=nothing, othercols=[], sortby=[], 
            rangenorm::Vector=[],
            percentnorm::Union{Nothing,Real}=nothing,
            unitnoskip::Bool=false,
            asmat::Bool=false)
        sfm = copy(sfm)
        stat = Symbol(stat)
        sortby = Symbol.(sortby)
        if filtval !== nothing
            if filtval === NaN
                sfm = sfm[(!).(isnan.(sfm[!,stat])), :]
            else
                sfm = sfm[sfm[!,stat] .!= filtval, :]
            end

        end
        rows = unique([:unit, othercols..., sortby...])
        if !(isempty(rangenorm))
            normfunc = x->DIutils.norm_extrema(x, rangenorm)
            sfm = combine(groupby(sfm, :unit), stat => normfunc => stat,
                         rows .=> rows, :shift => :shift)
        end
        if percentnorm !== nothing
            normfunc = x->DIutils.norm_percent(x, percentnorm)
            sfm = combine(groupby(sfm, :unit), stat => normfunc => stat,
                         rows .=> rows, :shift => :shift)
        end
        out = !(isempty(sortby)) ? sort(unstack(sfm, rows, :shift, stat), sortby) :
                                   unstack(sfm, rows, :shift, stat)
        if asmat
            unit = unitnoskip ? collect(1:length(out.unit)) : out.unit
            shifts = parse.(Float32, names(out[:, Not(rows)]))
            axs = [Axis{:unit}(unit), Axis{:shift}(shifts)]
            dat = Matrix(out[:, Not(rows)])
            #@info "ax" axs size(dat) shifts
            out = AxisArray(dat, axs...)
        else
            out
        end
    end

    """
        get_metmatrix

    Shorthand version of getstatmat
    """
	function get_metmatrix(sfsmet; met, srt, normrange)	
		tmpstatmat = shiftmetrics.getstatmat(sfsmet, met; 
                                             filtval=(met == "coherence" ? NaN : 0),
                                             asmat=true,
                                             unitnoskip=false,
                                             sortby=[srt],
                                             (normrange ? (;rangenorm=[0,1]) : (;percentnorm=0.5))...);
	end

end

