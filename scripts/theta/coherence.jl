using GoalFetchAnalysis
using MATLAB
mat"addpath(genpath('/usr/local/chronux_2_12/'))"
using Munge.lfp
using Serialization
using DrWatson

animal,day="RY16",36
lfp_ca1 = Load.load_lfp(animal,day,tet=:default)
lfp_pfc = Load.load_lfp(animal,day,tet="PFC")
coh     = coherence(lfp_ca1, lfp_pfc)

serialize(datadir("coherence.serial"), coh)
coh = deserialize(datadir("coherence.serial"))

