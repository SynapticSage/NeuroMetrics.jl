
function cellpath(animal::String, day::Int, tag::String=""; type="csv", kws...)
    if tag != "" && tag != "*"
        tag = "_$tag"
    end
    csvFile = DrWatson.datadir("exp_raw", "visualize_raw_neural",
                               "$(animal)_$(day)_cell$tag.$type")
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

function save_cells(cells::AbstractDataFrame, pos...; kws...)
    save_table(cells, pos...; tablepath=:cells, kws...)
end

function load_tetrode(animal, day)
    cells = load_cells(animal,day)
    groups = groupby(cells,"tetrode")
    tetrodes = DataFrame()
    for group = groups
        n_cells = size(group,1)
        row = DataFrame(group[1,:])
        row[!, :n_cells] .= n_cells;
        append!(tetrodes, row);
    end
    out = if "cell" in names(tetrodes)
            tetrodes[!, Not(:cell)]
        else
            tetrodes
        end
    return out
end

