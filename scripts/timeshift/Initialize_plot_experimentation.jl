using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))

# Grab our raw data
using DataFrames
using KernelDensity, Distributions
using Plots, Measures
using ProgressMeter
using StatsPlots
using DataFramesMeta
using DataStructures: OrderedDict
using Revise
using Distributions

using GoalFetchAnalysis
import Timeshift
import Field
import Utils
import Table
import Plot

operation     = Field.operation
model         = Field.model
recon_process = Field.recon_process
quickactivate(expanduser("~/Projects/goal-code"))
include(scriptsdir("fields","Include.jl"))
@assert Field.get_fields isa Function
@time spikes, beh, ripples, cells = Load.load("RY16", 36);

include(scriptsdir("timeshift", "TimeShift_setsOfInterest.jl"))

# --- GET TEH DATA ---
I = Timeshift.load_mains()
S = Timeshift.load_shuffles()
F = Timeshift.load_fields()
F = Utils.namedtup.lambda_values(F, x->x.Rₕ)
key = (;first(keys(I))..., datacut=:all)
#imax = sort(Timeshift.imax(Timeshift.info_to_dataframe(I[key], shift_scale=:minutes)))
#iall = sort(Timeshift.info_to_dataframe(I[key], shift_scale=:minutes))

# --- DATAFRAME-A-TIZE ---
@time I = Table.to_dataframe(I, key_name="shift",   explode=true)
@time S = Table.to_dataframe(S, key_name=:keyboard, explode=true)
@time F = Table.to_dataframe(F, key_name="shift",   explode=false)

# --- DESCRIBE ------
using Term: Panel
P = Panel( Panel("Basic stats   main↓  shuffle→") /
      Panel(string(describe(I)), width=100) *
      Panel(string(describe(S)), width=100)
   )
Q = Panel( Panel("Unique vars main↓  shuffle→") /
      Panel(string(describe(I[:,Not(:value)], unique=>:unique)), width=100) *
      Panel(string(describe(S[:,Not(:value)], unique=>:unique)), width=100)
   )
print(P / Q)
# ------


