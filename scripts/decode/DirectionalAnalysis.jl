using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
include(scriptsdir("decode","InitializeCairo.jl"))

(splitfig, split_num) = first(Iterators.product([true,false], 1:2))

thresh = Dict("likelihood"=>0.1, "acausal_posterior"=>0.99, "causal_posterior"=> 0.99)

dothresh  = false # DO THRESH IN LOADDATA.jl ... this should be false, I'm doing it below
dodisplay = false
spikecolorby = :time # :celltime | :time
sortby = :ranreltime # :celltime | :time
padding = [.15, .15]

@info split_num
global decode_file
global beh, spikes, cycles, dat, T, x, y
function get_color(x,cmap,nan_color,alpha_all=1)
    if isnan(x)
        nan_color
    else
        col=get(cgrad(cmap), x)
        RGBA(col.r, col.g, col.b, alpha_all)
    end
end


# Load data
#include(scriptsdir("decode","Initialize.jl"))
include(scriptsdir("decode","InitializeCairo.jl"))
decode_file=replace(decode_file, "split=0"=>"split=$split_num")
@info decode_file
@time include(scriptsdir("decode", "LoadData.jl"))
@time include(scriptsdir("decode", "PreprocessLFP.jl"))

ripples, cycles = annotate_vector_info(ripples, cycles, beh, lfp, dat, x, y, T)
lfp.phase_plot = utils.norm_extrema(lfp.phase, extrema(spikes.unit))
lfp.phase = utils.norm_extrema(lfp.phase, (-pi,pi))
lfp.raw      = Float32.(utils.norm_extrema(lfp.raw,      extrema(spikes.unit)))
lfp.broadraw = Float32.(utils.norm_extrema(lfp.broadraw, extrema(spikes.unit)))
lfp.cycle, lfp.phase = Int32.(lfp.cycle), Float32.(lfp.phase)
cycles = cycles[cycles.cycle.!=-1,:]

beh = annotate_pastFutureGoals(beh; doPrevPast=false)



