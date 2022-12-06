using GoalFetchAnalysis
using MATLAB
mat"addpath(genpath('/usr/local/chronux_2_12/'))"
using Munge.lfp
using Serialization
using DrWatson

animal,day="RY16",36

lfp_ca1 = Load.load_lfp(animal,day,tet=:default)
lfp_pfc = Load.load_lfp(animal,day,tet="PFC")
#coh     = coherence(lfp_ca1, lfp_pfc)
#serialize(datadir("coherence.serial"), coh)
coh = deserialize(datadir("coherence.serial"))

avgcoh = Munge.lfp._getavgcoh(;coh...)
dfa = Munge.lfp._getdfcoh(;avgcoh..., average=true)
#serialize(datadir("dfa.serial"),dfa)
#dfa = deserialize(datadir("dfa.serial"))
Load.save_table_at_path(dfa, datadir("exp_raw", "visualize_raw_neural", "RY16_36_avgcoh.arrow"), "arrow")
dfa = nothing
GC.gc()

df =  Munge.lfp._getdfcoh(;coh...)
#serialize(datadir("coherence.df.serial"), df)
Load.save_table_at_path(df, datadir("exp_raw", "visualize_raw_neural", "RY16_36_coh.arrow"), "arrow")
