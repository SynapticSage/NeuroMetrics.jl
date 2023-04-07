module cause

    using RecipesBase
    using Plots
    using Statistics, NaNStatistics
    using Infiltrator
    
    rect_default=true
    params_default = (;)
    setparams(p::NamedTuple) = @eval cause params_default = $p

    @userplot PlotMeanCause
    @recipe function plotmeancause(plt::PlotMeanCause; transform=identity, 
            timefact=1,
            rect=rect_default)
        if length(plt.args) == 1
            set_cause = plt.args[1]
        else
            time, set_cause = plt.args
        end
        set_cause = set_cause[[isassigned(set_cause,p) 
                               for p in eachindex(set_cause)]]
        set_cause = hcat(set_cause...)'
        m --> :circle
        fillrange --> 0
        label --> "thing a -> thing b"
        #@infiltrate
        set_cause = mean(transform(set_cause),dims=1)'
        set_cause = rect ? rectify(set_cause) : set_cause
        if length(plt.args) == 1
            time = Float64.(collect(1:length(set_cause))) .* timefact
        end
        time, set_cause
    end

    @userplot PlotMedianCause
    @recipe function plotmediancause(plt::PlotMedianCause; transform=identity, 
            timefact=1,
            rect=rect_default)
        if length(plt.args) == 1
            set_cause = plt.args[1]
        else
            time, set_cause = plt.args
        end
       set_cause = set_cause[[isassigned(set_cause,p) 
                               for p in eachindex(set_cause)]]
        set_cause = hcat(set_cause...)'
        m --> :circle
        fillrange --> 0
        label --> "thing a -> thing b"
        set_cause = median(transform(set_cause),dims=1)'
        set_cause = rect ? rectify(set_cause) : set_cause
        if length(plt.args) == 1
            time = Float64.(collect(1:length(set_cause))) .* timefact
        end
        time, set_cause
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
        if length(one_cause) < 30
            m --> :circle
        end
        fillrange --> 0
        label --> "thing a -> thing b"
        #@infiltrate
        rect ? rectify(transform(one_cause)) : transform(one_cause) 
    end

    function rectify(a::AbstractArray)
        @info "rectify"
        #a[a .< 0] .= 0
        a
    end

    # ============================
    # plot_alltimes.jl functions
    # ============================
    key_filter = nothing
    """
        setkeyfilter
    sets a module variable used to filter the results in a predictiveasmmetry
    dictionary. all of the downstream plotX functions use this filter key.

    if it's set to nothing (default), no filtration occurs
    """
    setkeyfilter(K::Union{AbstractDict,NamedTuple,Nothing}=nothing) = 
         @eval cause key_filter = $K
    """
        getcausedistovertime
    Designed to look into dict of predictive asymmetry results
    and return values by time
    Returns matrix of (Xᵢ, Yⱼᵢ) elements. 
        Xᵢ = timeᵢ = a sample index i scaled by 1/30 (our camera framerate)
        Yⱼᵢ is a value for value j and sample i
    """
    function getcausedistovertime(C::Dict)
        V = [V for V in values(C) if !ismissing(V)]
        [(i * 1/30, V[j][i]) for i in 1:last(params[:horizon]) for j in 1:length(V)]
    end
    """
        getmedian
    per behavior key Kᵦ ∈ set_cause::Dict, we have a dictionary of values for the different
    subsets of the data. This takes the *MEDIAN* of those subsets for each Kᵦ.
    """
    function getmedian(set_cause)
        set_cause = collect(values(set_cause))
        inds = [isassigned(set_cause,p) && !ismissing(p) for p in eachindex(set_cause)]
        set_cause = skipmissing(set_cause[inds])
        set_cause = hcat(set_cause...)'
        #vec(median(transform(set_cause),dims=2))
        vec(median((set_cause),dims=1))
    end
    """
        getmean
    per behavior key Kᵦ ∈ set_cause::Dict, we have a dictionary of values for the different
    subsets of the data. This takes the *MEAN* of those subsets for each Kᵦ.
    """
    function getmean(set_cause)
        set_cause = collect(values(set_cause))
        inds = [isassigned(set_cause,p) && !ismissing(p) for p in eachindex(set_cause)]
        set_cause = skipmissing(set_cause[inds])
        set_cause = hcat(set_cause...)'
        #vec(mean(transform(set_cause),dims=2))
        try
         isempty(set_cause) ? missing : vec(mean((set_cause),dims=1))
        catch; 
            @infiltrate; 
        end
    end
    """
        getdiff
    transforms the data at the bottom of the Dict tree to a diff instead of
    predictive asymmetry's cumulative form
    """
    function getdiff(X::AbstractDict)
        typeof(X)(k=>getdiff(v) for (k,v) in X)
    end
    function getdiff(X::AbstractArray)
        [diff(X); 0]
    end
    function getdiff(X::Missing)
        missing
    end
    """
        plotmedianplushist
    Standard plot, gives the MEDIAN curve for each Kᵦ, and to summarize the variation, gives a histogram
    of the samples that contribute to each median
    """
    function plotmedianplushist(C::Dict;ch=:black,cmc=:black,labelmc="",
    histalpha=0.6,params=params_default,kws...)
        if key_filter !== nothing
            C = Dict(k=>v for (k,v) in C if k ∈ key_filter)
        end 
        V = [V for V in values(C) if !ismissing(V)]
        histme = [(i * 1/30,V[j][i]) for i in 1:last(params[:horizon]) for j in 1:length(V)]
        histogram2d(histme;kws...,alpha=histalpha)
        plotmediancause!(V; timefact=1/30, c=cmc, markersize=1,label=labelmc,
                        fill=0, fillalpha=0.35, linewidth=0.2)
        hline!([0], c=ch, linewidth=3, linestyle=:dash, 
               xlabel="time (s)", ylabel="ℙ_asymmetry", label="")
    end
    function bin_the_curves(curves...; time=time, 
            bins=[(0,0.20), (0.20,1), (1,3.5)], statfun=mean)
        curves = collect(skipmissing(curves))
        inds = [Utils.in_range(time,bin) for bin in bins]
        avgs = zeros(length(curves), length(bins))
        for (c,curve) in enumerate(curves)
            for (b,ind) in enumerate(inds)
                avg = statfun(curve[ind])
                avgs[c,b] = avg
            end
        end
        avgs
    end
    # Jacknife summaries
    func_full = x->getmean(x)
    func_bin = x->mean(bin_the_curves.(x),dims=1)
    export leaveoneout
    """
        leaveoneout
    "Leave one out" function for jacknifing
    """
    function leaveoneout(D::AbstractDict; func)
        leaveoneout(collect(values(D)); func)
    end
    function leaveoneout(V::Array; func)
        V = collect(skipmissing(V))
        results = []
        for i in 1:length(V) 
            V′ = V[Not(i)]
            result = func(V′)
            result = if result isa Array && length(result) == 1
                result[1]
            else
                result
            end
            push!(results,result)
        end
        results
    end

end
