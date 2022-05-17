
function cellpath(animal::String, day::Int, tag::String=""; kws...)
    if tag != "" && tag != "*"
        tag = "_$tag"
    end
    csvFile = DrWatson.datadir("exp_raw", "visualize_raw_neural",
                               "$(animal)_$(day)_cell$tag.csv")
end
function load_cells(pos...; kws...)
    path = cellpath(pos...; kws...)
    if occursin("*", path)
        base, dir = basename(path), dirname(path)
        @debug "base=$base, dir=$dir"
        paths = glob(base, dir)
    else
        paths = [path]
    end

    cells = DataFrame()
    @showprogress 0.1 "loading cell files" for path in paths
        cell = CSV.read(path, DataFrame;
                 strict=false, missingstring=["NaN", "", "NaNNaNi"],
                 csvkws...)
        cells = isempty(cells) ? cell : outerjoin(cells, cell, on=:unit, makeunique=true)
        table.clean_duplicate_cols(cells)
    end
    return cells
end
function save_cells(cells::DataFrame, pos...; merge_if_exist::Bool=true, kws...)
    #kws = (;kws...) # satellite default is true
    csvFile = cellpath(pos...; kws...)
    #if merge_if_exist && isfile(csvFile)
    #    prevcells = load_cells(pos...; kws...)
    #    cells = outerjoin(cells, prevcells, makeunique=true)
    #    # TODO function to accept left or right dups
    #end
    println("Saving cell data at $csvFile")
    cells |> CSV.write(csvFile)
end


