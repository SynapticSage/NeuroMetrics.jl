using GoalFetchAnalysis
using DataFrames, DataFramesMeta
using ProgressMeter
using Infiltrator

animal,day = "RY16", 36
dfa = Load.load_avgcoh(animal,day)
Load.register(beh, dfa; transfer=["correct","cuemem"], on="time")
dfa[!,:C] = replace(dfa.C, NaN=>missing)
dropmissing!(dfa,:C)

function clean_cuecorr(df)
    df = dropmissing(df, [:cuemem, :correct])
    @subset(df, vec(any(:cuemem .== [0, 1]',dims=2)) .&& vec(any(:correct .== [0, 1]',dims=2)) .|| 
            (:cuemem .== -1 .&& :correct .== -1))
end
function get_effects(dfa,pos...)
    dfa = !isempty(pos) ? subset(dfa, pos...) : dfa
    groups = groupby(dfa, ["correct","cuemem"])
    effect  = dropmissing(sort(combine(groups, 
                                       :C=>nanmean, :S1=>nanmean, :S2=>nanmean, 
                                       nrow=>:rows, 
                                       :C=>(x->std(x)./sqrt(length(x)))=>:Cerr,
                                       :S1=>(x->std(x)./sqrt(length(x)))=>:S1err,
                                       :S2=>(x->std(x)./sqrt(length(x)))=>:S2err,
                                      ), 
                               [:cuemem,:correct]), :cuemem)
    clean_cuecorr(effect)
end

theta = get_effects(dfa, :freq=>f-> Utils.in_range(f, [6,12]))
delta = get_effects(dfa, :freq=>f-> Utils.in_range(f, [1, 3]))
ripple = get_effects(dfa, :freq=>f-> Utils.in_range(f, [150,220]))
nottheta = get_effects(dfa, :freq=>f-> Utils.not_in_range(f, [4,15]))
#not_in_theta_harmons = TODO
overall = get_effects(dfa)

@time spikes, beh, ripples, cells = Load.load(animal, day);
lfp_ca1 = Load.load_lfp(animal,day,tet=:default, subtract_earlytime=true)
lfp_ca1 = Munge.lfp.annotate_cycles(lfp_ca1)

spikes = Munge.spiking.isolated(spikes, lfp_ca1, refreshcyc=true)
dfa_theta  = combine(groupby(subset(dfa,:freq=>f->Utils.in_range(f,[6,12])), :time),
                     :C=>mean=>:Ctheta,renamecols=false)
Load.register(dfa_theta, spikes; on="time", transfer=["Ctheta"])
dfa_delta  = combine(groupby(subset(dfa,:freq=>f->Utils.in_range(f,[1,3])), :time),
                     :C=>mean=>:Cdelta,renamecols=false)
Load.register(dfa_delta, spikes; on="time", transfer=["Cdelta"])
dfa_ripple  = combine(groupby(subset(dfa,:freq=>f->Utils.in_range(f,[150,220])), :time),
                     :C=>mean=>:Cripple,renamecols=false)
Load.register(dfa_ripple, spikes; on="time", transfer=["Cripple"])
Load.register(beh, spikes; on="time", transfer=["cuemem","correct"])

# RESULT :: Isolated spikes more common in high theta coherence
iso = combine(groupby(spikes, :isolated), 
        :Ctheta=>c->nanmean(c[(!).(ismissing.(c))]),
        :Ctheta=>(c->(c=c[(!).(ismissing.(c))];nanstd(c)./length(c)))=>:Ctheta_err,
        :Cdelta=>c->nanmean(c[(!).(ismissing.(c))]),
        :Cdelta=>(c->(c=c[(!).(ismissing.(c))];nanstd(c)./length(c)))=>:Cdelta_err,
        :Cripple=>c->nanmean(c[(!).(ismissing.(c))]),
        :Cripple=>(c->(c=c[(!).(ismissing.(c))];nanstd(c)./length(c)))=>:Cripple_err,
        renamecols=false 
       )
clean_cuecorr(iso)

# CONTENTION :: Less of an effect when looked at split by task variables 
#                (task explains more than iso)
@time isocc = combine(groupby(spikes, [:isolated,:cuemem,:correct]), 
        :Ctheta=>c->nanmean(c[(!).(ismissing.(c))]),
        :Ctheta=>(c->(c=c[(!).(ismissing.(c))];nanstd(c)./length(c)))=>:Ctheta_err,
        :Cdelta=>c->nanmean(c[(!).(ismissing.(c))]),
        :Cdelta=>(c->(c=c[(!).(ismissing.(c))];nanstd(c)./length(c)))=>:Cdelta_err,
        :Cripple=>c->nanmean(c[(!).(ismissing.(c))]),
        :Cripple=>(c->(c=c[(!).(ismissing.(c))];nanstd(c)./length(c)))=>:Cripple_err,
        renamecols=false 
       )
clean_cuecorr(isocc)

sort!(dfa, [:time, :freq])

# Spike triggered average

get_matrix(df::DataFrame, col) = 
     Matrix(unstack(df[!,:], :time, :freq, col)[:,Not(:time)])
get_unstack(df::DataFrame, col) = unstack(df[!,:], :time, :freq, col)
unstack_dict = Dict(col=>get_unstack(dfa, col) for col in [:C,:S1,:S2,:phi])

function spike_triggered_average(spikes::DataFrame, dfa::DataFrame)
    inds = Utils.searchsortednearest.([dfa.time], spikes.time)
    inds = inds[ abs.(dfa.time[inds] .- spikes.time) .< 0.05 ]
    times = dfa.time[inds]
    #dfg  = groupby(dfa,:time)
    #kT   = keys(dfg)
    #kTv  = [v[1] for v in kT]
    #freqs = unique(dfa.freq)
    #nF = length(freqs)

    xvals,yvals = Vector{Union{Missing,Vector}}(missing,length(times)),
                 Vector{Union{Missing,Matrix}}(missing,length(times))
    prog = Progress(length(inds);desc="triggering")
    Threads.@threads for (i,time) in collect(enumerate(times))
        
        #@infiltrate

        #keyloc = findfirst( time .== kTv )
        #I = UnitRange(clamp.(100 .* (-1,1) .+ keyloc, [1], [size(dfg,1)])...)
        #M = get_unstack(combine(dfg[kT[I]],identity), prop)
        
        df = subset(dfa, :time=>t->Utils.in_range(t, time .+ [-1, 1]))

        xvals[i] = df.time
        yvals[i] =  Matrix(df[!,Not(:time)])
        next!(prog)
    end
    (;times=xvals, yvals)
end

# Spectral averaging
using Serialization
Cta = Dict()
freq=unique(dfa.freq);

for (iso, area, field) in Iterators.product((true,false),("CA1","PFC"),(:C,))
    key = (;iso, area, field)
    if key ∉ keys(Cta)
        Cta[key] = spike_triggered_average(
                                @subset(spikes,:isolated .== iso .&& :area .== area),
                                        unstack_dict[field])
        #serialize(datadir("checkpoint"),(;Cta, freq))
    end
end

serialize(datadir("checkpoint"),(;Cta,freq))

using Statistics, NaNStatistics
using Plots
function get_plot(key)
    cta = Cta[key]
    sel = length.(cta.times) .== median(length.(cta.times))
    times = mean(cta.times[sel])
    times = times .- minimum(times)
    vals  = mean(cta.yvals[sel])
    heatmap(times,freq,vals',
            c=:vik, size=1/2 .*(800,2000))
end
h1 = get_plot((iso=true,  area = "CA1", field = :C))
h2 = get_plot((iso=true,  area = "PFC", field = :C))
h3 = get_plot((iso=false, area = "CA1", field = :C))
h4 = get_plot((iso=false, area = "PFC", field = :C))
plot(h1,h2, h3, h4, size=(1600,2000))
