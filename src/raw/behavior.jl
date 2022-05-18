
function behaviorpath(animal::String, day::Int, tag::String="")
    tag  = length(tag) == 0 ? tag : "_$tag"
    path = datadir("exp_raw", "visualize_raw_neural",
                     "$(animal)_$(day)_beh$tag")
    @debug "path=$path"
    if occursin("*",tag)
        path = glob(basename(path), dirname(path))
    end
    @debug "path=$path"
    return path
end

function load_behavior(animal::String, day::Int, tag::String="")
    function typeFunc(type, name)
        if occursin("Vec", string(name))
            type = ComplexF32;
        elseif name == "time" type = Float32;
        else
            type = nothing;
        end
        return type
    end
    typemap = Dict(Int64=>Int16);
    @debug "animal=$animal, day=$day, tag=$tag"
    behCSV = behaviorpath(animal, day, tag) .* ".csv"
    @info behCSV
    readFile(file) = CSV.read(file, DataFrame;
                              strict=false, 
                              missingstring=["NaNNaNi", "NaNNaNi,", ""], 
                              types=typeFunc, typemap=typemap, csvkws...)
    if behCSV isa Vector
        beh = [readFile(file) for file in behCSV]
        beh = hcat(beh...)
    else
        beh = readFile(behCSV)
    end
    if beh.time isa Vector{String}
        beh.time = parse.(Float32, beh.time);
    end
    for col in names(beh)
        if occursin("current", col) || occursin("egoVec", col)
            replace!(beh[!,col], missing=>NaN)
        end
    end
    @assert ("x" ∈ names(beh)) "Fuck"
    return beh
end

function save_behavior(behavior::AbstractDataFrame, pos...; kws...)
    save_table(behavior, pos...; tablepath=:behavior, kws...)
end

module behavior
    function add_next_target(beh)
    end
    function add_previous_target(beh)
    end
    additions = Dict("previous_target" => add_previous_target,
                     "next_target"=> add_next_target)
end
export behavior

