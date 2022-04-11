quickactivate("/home/ryoung/Projects/goal-code/")
includet(srcdir("field.jl"))
includet(srcdir("filt.jl"))
includet(srcdir("raw.jl"))
includet(srcdir("field/timeshift.jl"))
includet(srcdir("field/info.jl"))
spikes, beh = raw.load("RY16", 36, data_source=["spikes","behavior"])
props = ["x", "y"]
splitby=["unit", "area"]
kws=(;splitby, filters=merge(filt.speed_lib, filt.cellcount))
newkws = (; kws..., resolution=40, gaussian=2.3*0.5, props=props,
          filters=merge(kws.filters))

# SINGLE TIME SHIFT
# 6-20 seconds for 1
@time place = timeshift.get_field_shift(beh, spikes, 1; newkws...);

# MULTIPLE TIME SHIFT
#
# -----------
# 16 threads
# -----------
# (run 1) 400 seconds for 20, with multi-threading, with 90% of that compile time
# (run 2) 180 seconds 85% compile time
# (run 3) 131 seconds 76% compile time
# (run 4) 131 seconds 76% compile time
#
# ---------
# 4 threads
# ---------
#
@time place = timeshift.get_field_shift(beh, spikes, -1:0.1:1; multithread=true, newkws...);

# ---------
# 1 thread
# ---------
# (run 1) 108 seconds for 20 shifts without multi-threading, with 10% compile time
@time place = timeshift.get_field_shift(beh, spikes, -1:0.1:1; multithread=false, newkws...);
F["place"] = place
