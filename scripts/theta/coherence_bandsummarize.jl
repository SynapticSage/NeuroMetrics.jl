using Statistics, NaNStatistics, Bootstrap
using StatsPlots
using GoalFetchAnalysis
using DataFrames, DataFramesMeta

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
