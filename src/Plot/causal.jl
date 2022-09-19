module cause

    using RecipesBase
    using Plots
    using Statistics, NaNStatistics
    
    rect_default=false

    @userplot PlotMeanCause
    @recipe function plotmeancause(plt::PlotMeanCause; transform=identity, rect=rect_default)
        set_cause = plt.args[1]
        set_cause = set_cause[[isassigned(set_cause,p) 
                               for p in eachindex(set_cause)]]
        set_cause = hcat(set_cause...)'
        m --> :circle
        fillrange --> 0
        label --> "thing a -> thing b"
        #@infiltrate
        set_cause = mean(transform(set_cause),dims=1)'
        set_cause = rect ? rectify(set_cause) : set_cause
        set_cause
    end

    @userplot PlotMeanCauseDiff
    @recipe function plotmeancause(plt::PlotMeanCauseDiff; transform=identity, rect=rect_default)
        set_cause = plt.args[1]
        set_cause = set_cause[[isassigned(set_cause,p) 
                               for p in eachindex(set_cause)]]
        set_cause = [diff(p) for p in set_cause]
        set_cause = hcat(set_cause...)'
        m --> :circle
        fillrange --> 0
        label --> "thing a -> thing b"
        #@infiltrate
        set_cause = mean(transform(set_cause),dims=1)'
        set_cause = rect ? rectify(set_cause) : set_cause
        set_cause
    end

    @userplot PlotCauseDiff
    @recipe function plotcausediff(plt::PlotCauseDiff; transform=identity, rect=rect_default)
        one_cause = plt.args[1]
        diff_one_cause = diff(one_cause)
        m --> :circle
        fillrange --> 0
        label --> "thing a -> thing b"
        #@infiltrate
        rect ? transform(diff_one_cause) : rectify(transform(diff_one_cause))
    end

    @userplot PlotCause
    @recipe function plotcause(plt::PlotCause; transform=identity, rect=rect_default)
        one_cause = plt.args[1]
        m --> :circle
        fillrange --> 0
        label --> "thing a -> thing b"
        #@infiltrate
        rect ? transform(one_cause) : rectify(transform(one_cause))
    end

    function rectify(a::AbstractArray)
        a[a .< 0] .= 0
        a
    end

end
