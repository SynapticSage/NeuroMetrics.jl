module Decode

    using DrWatson
    using Reexport

    using Statistics
    using LoopVectorization
    using Revise
    __revise_mode__ = :evalassign

    using DIutils
    import Load

    searchsortednearest = Utils.searchsortednearest
    export searchsortednearest

    function movingmean(dat)
        R = [rolling(mean, dat[i,j,:],5) for i in 1:size(dat,1), j in 1:size(dat,2)]
        datnew = zeros(size(dat,1),size(dat,2), length(R[1]))
        for i in 1:size(dat,1), j in 1:size(dat,2)
            datnew[i,j,:] = vec(R[i,j])
        end
        dat=datnew
    end

    function quantile_threshold(dat, thresh=nothing; sample_dim=3,
                                                     nan_replace_val=1)
        m = [quantile(Utils.skipnan(vec(d)), thresh) 
             for d in eachslice(dat, dims=sample_dim)]
        m = reshape(m, (1,1,length(m)))
        sz = size(dat)
        dat[dat .<= m] .= NaN
        return reshape(dat, sz)
    end

    include(srcdir("Decode","checkpoint.jl"))
    @reexport using .checkpoint
    include(srcdir("Decode","makie_observable.jl"))
    @reexport using .makie_observable
    include(srcdir("Decode","plot.jl"))
    @reexport import .plot
    include(srcdir("Decode","bayes.jl"))
    @reexport import .bayes


end
