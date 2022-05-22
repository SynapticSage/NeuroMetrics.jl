module decode
using Statistics
using LoopVectorization
using Revise
__revise_mode__ = :evalassign

include("utils.jl")
using .utils
include("raw.jl")
import .raw
include("table.jl")
import .table

searchsortednearest = utils.searchsortednearest
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
    m = [quantile(utils.skipnan(vec(d)), thresh) 
         for d in eachslice(dat, dims=sample_dim)]
    m = reshape(m, (1,1,length(m)))
    sz = size(dat)
    dat[dat .<= m] .= NaN
    return reshape(dat, sz)
end

include("decode/preprocess_lfp.jl")
include("decode/preprocess_beh.jl")
include("decode/checkpoint.jl")
include("decode/makie_observable.jl")
include("decode/plot.jl")
#include("decode/wellanalysis.jl")

end
