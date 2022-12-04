using GoalFetchAnalysis
using Munge.lfp
animal,day="RY16",36
lfp_ca1 = Load.load_lfp(animal,day,tet=:default)
lfp_pfc = Load.load_lfp(animal,day,tet="PFC")
coh = coherence(lfp_ca1, lfp_pfc)
