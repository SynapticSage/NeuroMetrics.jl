
function taskpath(animal, day)
    csvFile = DrWatson.datadir("exp_raw",
                               "visualize_raw_neural",
                               "$(animal)_$(day)_task.csv")
end

function load_task(animal::String, day::Int)
    taskFile = taskpath(animal, day)
    typemap = Dict(Int64=>Int16);
    @info csvFile
    task = CSV.read(csvFile, DataFrame; strict=false, typemap=typemap,
             missingstring=["NaN", "NaNNaNi", "NaNNaNi,", ""], csvkws...)
    transform!(task, :x => pxtocm => :x, :y => pxtocm => :y)
    return task
end

function save_task(behavior::AbstractDataFrame, pos...; kws...)
    save_table(behavior, pos...; tablepath=:task, kws...)
end

function well_locations(task; day=nothing, epoch=nothing)
end
