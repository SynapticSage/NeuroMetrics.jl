using Statistics, NaNStatistics, Bootstrap
using StatsPlots
using GoalFetchAnalysis
using DataFrames, DataFramesMeta

animal,day = "RY16", 36
dfa = Load.load_avgcoh(animal,day)
Load.register(beh, dfa; transfer=["correct","cuemem"], on="time")
dfa[!,:C] = replace(dfa.C, NaN=>missing)
dropmissing!(dfa,:C)

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
    effect = @subset(effect, vec(any(:cuemem .== [0, 1]',dims=2)) .&& vec(any(:correct .== [0, 1]',dims=2)) .|| 
            (:cuemem .== -1 .&& :correct .== -1))
end

theta = get_effects(dfa, :freq=>f-> Utils.in_range(f, [6,12]))
delta = get_effects(dfa, :freq=>f-> Utils.in_range(f, [1, 3]))
ripple = get_effects(dfa, :freq=>f-> Utils.in_range(f, [150,220]))
nottheta = get_effects(dfa, :freq=>f-> Utils.not_in_range(f, [4,15]))
#not_in_theta_harmons = TODO
overall = get_effects(dfa)


