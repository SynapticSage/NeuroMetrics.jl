
function cellpath(animal::String, day::Int, tag::String=""; type="csv", kws...)
    if tag != "" && tag != "*"
        tag = "_$tag"
    end
    csvFile = DrWatson.datadir("exp_raw", "visualize_raw_neural",
                               "$(animal)_$(day)_cell$tag.$type")
end

function cellpaths(animal::String, day::Int, tag::String=""; type="csv", kws...)
    path = cellpath(animal, day, tag; kws...)
    if occursin("*", path)
        base, dir = basename(path), dirname(path)
        @debug "base=$base, dir=$dir"
        paths = glob(base, dir)
    else
        paths = [path]
    end
    return paths
end


function load_cells(pos...; kws...)
    paths  = cellpaths(pos...; kws...)
    cells = DataFrame()
    @showprogress 0.1 "loading cell files" for path in paths
        cell = CSV.read(path, DataFrame; csvkws...)
        cells = isempty(cells) ? cell : outerjoin(cells, cell, on=:unit, makeunique=true)
        table.clean_duplicate_cols(cells)
    end
    return cells
end

function save_cells(cells::AbstractDataFrame, pos...; kws...)
    save_table(cells, pos...; tablepath=:cells, kws...)
end
function save_cell_taginfo(cells::AbstractDataFrame, animal::String, day::Int, tag::String; kws...)
    save_table(cells, animal, day, tag; tablepath=:cells, kws...)
end

function cells_to_type(animal::String, day::Int, tag::String="*", 
        from::String="csv", to::String="arrow")
    paths = cellpaths(animal, day, tag)
    for path in paths
        data = load_table_at_path(path, from; load_kws...)
        path = replace(path, "."*from=>"."*to)
        savekws=(;)
        save_table_at_path(data, path, type; save_kws...)
    end

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
