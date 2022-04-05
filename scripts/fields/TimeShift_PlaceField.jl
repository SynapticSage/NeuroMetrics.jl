includet(srcdir("timeshift.jl"))
props = ["x", "y"]
newkws = (; kws..., resolution=resolution, gaussian=2.3*0.5, props=props,
          filters=merge(kws.filters))
place = timeshift.get_field_shift(beh, spikes, 1; newkws...);
place = timeshift.get_field_shift(beh, spikes, -1:0.1:1; newkws...);
F["place"] = place

