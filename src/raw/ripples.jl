
function load_ripples(animal, day)
    typemap = Dict(Int64=>Int16);
    rippleFile = DrWatson.datadir("exp_raw",
                               "visualize_raw_neural",
                               "$(animal)_$(day)_ripple.csv")
    @info rippleFile
    ripples = CSV.read(rippleFile, DataFrame; strict=false,
             typemap=typemap,
             missingstring=["NaN", "NaNNaNi", "NaNNaNi,", ""],
             csvkws...)
end

