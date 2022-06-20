@time include(scriptsdir("timeshift", "Initialize.jl"))
include(scriptsdir("timeshift", "TimeShift_setsOfInterest.jl"))

# --- GET TEH DATA ---
I = Timeshift.load_mains()
S = Timeshift.load_shuffles()
F = Timeshift.load_fields()
key = (;first(keys(I))..., datacut=:all)
#imax = sort(Timeshift.imax(Timeshift.info_to_dataframe(I[key], shift_scale=:minutes)))
#iall = sort(Timeshift.info_to_dataframe(I[key], shift_scale=:minutes))


# --- DATAFRAME-A-TIZE ---
@time I = Table.to_dataframe(I, key_name="shift", explode=true)
@time S = Table.to_dataframe(S, key_name=:keyboard, explode=true)
@time F = Table.to_dataframe(S, key_name="shift", explode=false)

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


