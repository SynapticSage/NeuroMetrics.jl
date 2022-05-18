
function ripplespath(animal, day)
    rippleFile = DrWatson.datadir("exp_raw",
                               "visualize_raw_neural",
                               "$(animal)_$(day)_ripple.csv") 
end

function load_ripples(animal, day)
    rippleFile = ripplespath(animal,day)
    typemap = Dict(Int64=>Int16);
    @info rippleFile
    ripples = CSV.read(rippleFile, DataFrame; strict=false,
             typemap=typemap,
             missingstring=["NaN", "NaNNaNi", "NaNNaNi,", ""],
             csvkws...)
end

function save_ripples(ripples::AbstractDataFrame, pos...; kws...)
    save_table(ripples, pos...; tablepath=:ripples, kws...)
end
