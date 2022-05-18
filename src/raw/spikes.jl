
function spikespath(animal::String, dat::Int)
    rawSpikingCSV = DrWatson.datadir("exp_raw",
                                     "visualize_raw_neural",
                                     "$(animal)_$(day)_labeled_spiking.csv"
                                    )
end

function load_spikes(animal::String, day::Int; beh=Nothing)
    @info rawSpikingCSV
    raster = CSV.read(rawSpikingCSV, DataFrame;
             strict=false, missingstring=["NaN", "", "NaNNaNi"],
             csvkws...)

    # And let's add some brain area specific numbering
    groups = groupby(raster, "area");
    for g = 1:length(groups)
        unit = unique(groups[g].unit)
        areaunit = 1:length(unit);
        groups[g].areaunit = map(x -> Dict(unit .=> areaunit)[x],
                                 groups[g].unit)
    end
    raster = combine(groups, x->x)
end

function save_behavior(spikes::AbstractDataFrame, pos...; kws...)
    save_table(spikes, pos...; tablepath=:spikes, kws...)
end
