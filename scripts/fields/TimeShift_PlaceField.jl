quickactivate("/home/ryoung/Projects/goal-code/")
includet(srcdir("field.jl"))
includet(srcdir("filt.jl"))
includet(srcdir("raw.jl"))
includet(srcdir("field/timeshift.jl"))
spikes, beh = raw.load("RY16", 36, data_source=["spikes","behavior"])
props = ["x", "y"]
splitby=["unit", "area"]
kws=(;splitby, filters=merge(filt.speed_lib, filt.cellcount))
newkws = (; kws..., resolution=40, gaussian=2.3*0.5, props=props,
          filters=merge(kws.filters))
@time place = timeshift.get_field_shift(beh, spikes, 1; newkws...);
@time place = timeshift.get_field_shift(beh, spikes, -1:0.1:1; newkws...);
F["place"] = place

