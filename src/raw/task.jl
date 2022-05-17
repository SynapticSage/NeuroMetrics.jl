
function load_task(animal::String, day::Int)
    typemap = Dict(Int64=>Int16);
    csvFile = DrWatson.datadir("exp_raw",
                               "visualize_raw_neural",
                               "$(animal)_$(day)_task.csv")
    @info csvFile
    task = CSV.read(csvFile, DataFrame; strict=false, typemap=typemap,
             missingstring=["NaN", "NaNNaNi", "NaNNaNi,", ""], csvkws...)
    transform!(task, :x => pxtocm => :x, :y => pxtocm => :y)
    return task
end

function well_locations(task; day=nothing, epoch=nothing)
end
