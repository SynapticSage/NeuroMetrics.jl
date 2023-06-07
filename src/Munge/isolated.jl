__precompile__(false)
"""
    isolated

This module contains code used to analyze isolated spikes. It is
designed to be used in conjunction with the `GoalFetchAnalysis` module
and the `GoalFetchAnalysis.Munge` module.
"""
module isolated


# General purpose packages
using DataFrames, DataFramesMeta, DrWatson, JLD2, Missings, ProgressMeter, Random, Statistics, ThreadSafeDicts

# Data structure and utility packages
using DIutils.binning
using DIutils.dict
using DataStructures: OrderedDict
using StatsBase

# Math and statistical packages
using Distributions, GLM, GLMNet, Infiltrator, Metrics, MultivariateStats, NaNStatistics, RecipesBase

# Plotting and visualization packages
using Plots

# Machine Learning packages
using MLJ, ScikitLearn

# Interfacing with other languages
using MATLAB
using MATLAB: matlab_cmd
using PyCall, RCall

# ArgParse package
using ArgParse
using ArgParse: test_range

# Specific imports
import DI, DIutils, Distributions, Logging, MATLAB
import DIutils: Table

# function __init__()
#     @eval isolated pyglmnet = pyimport("pyglmnet")
#     @eval isolated using MLJScikitLearnInterface: 
#          ElasticNetCVRegressor, LogisticCVClassifier
# end
                                                            
#  . .     ,---.|              |              o     |         
# -+-+-    |    |---.,---.,---.|__/ ,---.,---..,---.|--- ,---.
# -+-+-    |    |   ||---'|    |  \ |   ||   |||   ||    `---.
#  ` `     `---'`   '`---'`---'`   `|---'`---'``   '`---'`---'
# SECTION: NEW

"""
    import_dataframes(animal, day;
        checkpoint, lfp_convert, super_convert, BEH, SPIKES, animals)

#Arguments
- `animal`: string, name of animal
- `day`: integer, day of animal
- `checkpoint`: boolean, whether to load from checkpoint
- animals: vector of tuples, each tuple is (animal, day)
# Returns
- `cells`: DataFrame, cell data
- `spikes`: DataFrame, spike data
- `beh`: DataFrame, behavior data, if available, otherwise nothing
- `ripples`: DataFrame, ripple data
- `cycles`: DataFrame, cycle data
- `lfp`: DataFrame, lfp data, if available
"""
function import_dataframes(animal::String, day::Int;
    tetrode_set::Union{Symbol,String}="best",
    checkpoint::Bool=true, lfp_convert=nothing, super_convert=nothing)
    lfp_convert   = lfp_convert === nothing ? DI.get_0time_pre_superanimal() : lfp_convert
    super_convert = super_convert === nothing ? DI.convert_super_to_time0() : super_convert
    # ONE ANIMAL CHECKPOINT
    if checkpoint && !occursin("super", animal)
        println("Loading from checkpoint for $animal, day $day")
        lfp, cycles, spikes, ripples, beh = get_animal_checkpoint(animal, day; tetrode_set)
        cells = DI.load_cells("super_clean", 0)
        cells = subset(cells, :animal => a->a.==animal, :day => d->d.==day)
    # MULTIPLE ANIMALS CHECKPOINT
    elseif checkpoint
        println("Loading from checkpoint for $(DI.animaldays())")
        OUT = []
        animals, days = DI.animaldays()
        for (animal, day) in zip(animals, days)
            lfp, cycles, spikes, ripples, beh = get_animal_checkpoint(animal, day; tetrode_set)
            push!(OUT, (lfp, cycles, spikes, ripples, beh))
        end
        lfp     = vcat([x[1] for x in OUT]...)
        cycles  = vcat([x[2] for x in OUT]...)
        spikes  = vcat([x[3] for x in OUT]...)
        ripples = vcat([x[4] for x in OUT]...)
        beh  = vcat([x[5] for x in OUT]...)
        cells = DI.load_cells("super_clean", 0)
    # ONE ANIMAL from scratch
    else
        println("Loading from scratch for $animal, day $day")

        BEH, SPIKES, cells = DI.load("super_clean", 0; data_sources=["beh", "spikes", "cells"])
        println("lfp_convert: ", lfp_convert)
        println("super_convert: ", super_convert)
        beh = subset(BEH, :animal => a->a.==animal, :day => d->d.==day)
        animals = unique(beh.animal)
        for animal in animals
            println("Correcting $animal by ", super_convert[animal])
            BEH[BEH.animal         .== animal, :time] .-= super_convert[animal]
            SPIKES[SPIKES.animal   .== animal, :time] .-= super_convert[animal]
        end
        spikes = subset(SPIKES, :animal => a->a.==animal, :day => d->d.==day, view=true)
        ripples = DI.load_ripples(animal, day)
        ripples = DataFrames.transform(ripples,
                    :time => t-> t .- lfp_convert[animal],
                    :start => t-> t .- lfp_convert[animal],
                    :stop => t-> t .- lfp_convert[animal],
                    renamecols=false)
        lfp = cycles = nothing
    end

    return cells, spikes, beh, ripples, cycles, lfp
end


"""
    get_animal_checkpoint(animal, day; tetrode_set="best")

Returns a tuple of DataFrames containing the checkpointed data
used for iso analyses.

# Arguments
- `animal::String`: The animal name
- `day::Int`: The day number
- `tetrode_set::String`: The tetrode set to use. Defaults to `"best"`.
"""
function get_animal_checkpoint(animal::String, day::Int; tetrode_set::String="best")
    println("Loading $tetrode_set data for $animal")
    zero_of_beh = DI.get_0time_pre_superanimal()
    lfp = DI.load_lfp(animal, day; append="$tetrode_set")
    lfp = DataFrames.transform(lfp, :time=> t-> t .- zero_of_beh[animal],
        renamecols=false)
    # DataFrames.transform!(lfp, :time=> t-> t .+ time_factors[animal], renamecols=false)
    cycles       = DI.load_cycles(animal, day, "$tetrode_set")
    cycles       = DataFrames.transform(cycles,
        :start=> t-> t .- zero_of_beh[animal],
        :stop=> t-> t .-  zero_of_beh[animal],
        renamecols=false
    )
    spikes = DataFrames.transform(DI.load_spikes(animal, day, "$(tetrode_set)_cycles_isolated"),
        :time=> t-> t .- zero_of_beh[animal],
        renamecols=false
    )
    # spikes = subset(spikes, :animal=>a->a.==animal, :day=>d->d.==day)
    # Print out extrema(time) for each of these just to be sure
    convert_to_f32 = [:broadraw  :phase :ripple  :rippleamp  :ripplephase]
    for col in convert_to_f32
        if hasproperty(lfp, col)
            lfp[!,col] = convert(Array{Float32,1}, lfp[!,col])
        end
    end
    ripples = DI.load_ripples(animal, day)
    ripples = DataFrames.transform(ripples, 
        :time=> t-> t .- zero_of_beh[animal],
        :start=> t-> t .- zero_of_beh[animal],
        :stop=> t-> t .- zero_of_beh[animal],
        renamecols=false)
    beh = DataFrames.transform(DI.load_behavior(animal, day),
        :time=> t-> t .- zero_of_beh[animal],
        renamecols=false)
    return lfp, cycles, spikes, ripples, beh
end

#  . .     ,---.|              |              o     |         
# -+-+-    |    |---.,---.,---.|__/ ,---.,---..,---.|--- ,---.
# -+-+-    |    |   ||---'|    |  \ |   ||   |||   ||    `---.
#  ` `     `---'`   '`---'`---'`   `|---'`---'``   '`---'`---'
#                                   |                         
# SECTION: OLD

function load_checkpoint()
end

export path_iso
"""
    path_iso(animal::String, day::Int, tet=:ca1ref)::String

Returns the path to the isolated data for the given animal, day, and
tet. The default tet is the ca1ref tet.
"""
function path_iso(animal::String, day::Int, tet=:ca1ref)::String
    datadir("isolated","iso_animal=$(animal)_day=$(day)_tet=$(tet).jld2")
end
"""
    path_iso(opt::AbstractDict)::String

Returns the path to the isolated data for the given animal, day, and
tet. The default tet is the ca1ref tet.
"""
function path_iso(opt::AbstractDict)::String
    path_iso(opt["animal"], opt["day"], opt["tet"])
end
"""
    path_iso(pos...; append::String="_cyclewise")::String
Returns the path to the isolated data for the given animal, day, and
tet. The default tet is the ca1ref tet. Optionally, append a string.
"""
function path_iso(pos...; append::String)::String
    f = path_iso(pos...)
    replace(f, ".jld2" => "$(append).jld2")
end

export load_iso
"""
    load_iso(pos...)

# Arguments
----------
see also: path_iso for possible pos arguments

# Returns
-------
OrderedDict
    keys are the variable names and values are the variables

see also: path_iso
"""
function load_iso(pos...)::OrderedDict
    results = OrderedDict()
    storage=JLD2.jldopen(path_iso(pos...),"r")
    try
        results = OrderedDict(zip(keys(storage),
                                  [storage[k] for k in keys(storage)]))
    catch exception
        throw(exception)
    finally
        close(storage)
    end
    results
end
"""
    load_iso

Modifies the module to incldude isolated variables loaded

see also: path_iso for possible pos arguments
"""
function load_iso!(Mod::Module, pos...)::Nothing
    data = load_iso(pos...)
    dict.load_keysvals_to_module!(Mod, keys(data), values(data))
end

export parser
"""
    parse(args=nothing; return_parser::Bool=false)

return command line parser for flags that control my manifold
related analyses
"""
function parser(args=nothing; return_parser::Bool=false)
    parser = ArgParseSettings()
    @add_arg_table parser begin
        "--animal"
            help = "the animal to run default: the super animal"
            arg_type = String
            default = "RY16"
        "--day"
            help = "the day to run"
            arg_type = Int
            default = 36
        "--tet"
            help = "tetrode to use for isolated spikes"
            arg_type = String
            default = "ca1ref"
        "--dataset"
             help = "dataset preset"
             arg_type = Int
             default = 0
        "--filt", "-f"
            help = "inactive for now"
            arg_type = Symbol
            default = :all
        "--cycles", "-c"
            help = "how many cycles to explore ahead/behind"
            arg_type = Int
            default = 8
        "--overwrite"
            help = "overwite checkpointed data"
            arg_type = Bool
            action = "store_true"
            default = false
    end
    if return_parser
        return parser
    else
        if args !== nothing
            opt = postprocess(parse_args(args, parser))
        else
            opt = postprocess(parse_args(parser))
        end
        @info "parsed options" opt
        return opt
    end
end

function postprocess(opt::AbstractDict)
    if tryparse(Int64, opt["tet"]) isa Number
        opt["tet"] = Base.parse(Int64, opt["tet"])
    else
        opt["tet"] = Symbol(opt["tet"])
    end
    opt
end

export construct_predict_isospikecount
"""
    construct_predict_isospikecount(df, cells, input_area="CA1";
            other_vars=[], other_ind_vars=[])

Returns a list of GLM formulae for predicting the spike count of
each neuron in the dependent area from the spike count of each
neuron in the independent area. The independent area is assumed
to contain isolated spike columns.

# Inputs
- `df`: the dataframe containing the spike counts
- `cells`: the dataframe containing the cell information
- `input_area`: the area that contains the isolated spikes
- `other_vars`: other variables to include in the model
- `other_ind_vars`: other variables to include in the model as
    independent variables

# Returns
- `formulae`: a list of StatsModels formulae
"""
function construct_predict_isospikecount(df, cells, ind_area="CA1";
        dep_area=nothing, other_vars=[], other_ind_vars=[])
    uArea = unique(cells.area)
    @assert length(uArea) == 2 "Only supports two area dataframes"
    dep_area = dep_area === nothing ? 
        dep_area = setdiff(uArea, [ind_area]) : dep_area
    dep_neurons = @subset(cells,:area .==dep_area).unit
    ind_neurons = @subset(cells,:area .==ind_area).unit
    dep_neurons = string.(dep_neurons) .* "_i"
    ind_neurons = string.(ind_neurons) 
    filter!(n->n ∈ names(df), dep_neurons)
    filter!(n->n ∈ names(df), ind_neurons)
    ind_eq_dep = ind_area == dep_area
    
    formulae = []
    for n_dep in dep_neurons
        ind_neuron_set = ind_eq_dep ? 
            setdiff(ind_neurons, [n_dep]) : ind_neurons 
        n_ind = first(ind_neurons)
        independents = GLM.Term(Symbol(n_ind))
        for n_ind in ind_neuron_set[2:end]
            independents += GLM.Term(Symbol(n_ind)) 
        end
        formula = GLM.Term(Symbol(n_dep)) ~ independents
        push!(formulae, formula)
    end
    Vector{FormulaTerm}(formulae)
end

export construct_predict_spikecount
"""
    construct_predict_spikecount(df, cells, indep_area="CA1";
            dep_area=nothing, other_vars=[], other_ind_vars=[])

Returns a list of GLM formulae for predicting the spike count of
each neuron in the dependent area from the spike count of each
neuron in the independent area. The independent area is raw spike
counts.

# Inputs
- `df`: the dataframe containing the spike counts
- `cells`: the dataframe containing the cell information
- `indep_area`: the area that contains spikecount to predict from
- `dep_area`: the area that contains spikecount to predict

# Returns
- `formulae`: a list of StatsModels formulae
"""
function construct_predict_spikecount(df, cells::DataFrame, 
        indep_area="CA1"; dep_area=nothing, 
        other_vars=[], other_ind_vars=[])

    uArea = unique(cells.area)
    @assert length(uArea) == 2 "Only supports two area dataframes"
    dep_area = dep_area === nothing ? 
                setdiff(uArea, [indep_area]) : dep_area
    ind_eq_dep = indep_area == dep_area

    dep_neurons = @subset(cells,:area .== dep_area).unit
    ind_neurons = @subset(cells,:area .== indep_area).unit
    filter!(n->string(n) ∈ names(df), dep_neurons)
    filter!(n->string(n) ∈ names(df), ind_neurons)
    
    formulae = Vector{FormulaTerm}()
    for n_dep in dep_neurons
        ind_neuron_set = ind_eq_dep ? 
            setdiff(ind_neurons, [n_dep]) : ind_neurons 
        n_ind = first(ind_neuron_set)
        independents = GLM.Term(Symbol(n_ind))
        for ni in ind_neuron_set[2:end]
            independents += GLM.Term(Symbol(ni)) 
        end
        formula = GLM.Term(Symbol(n_dep)) ~ independents
        push!(formulae, formula)
    end
    formulae
end

export construct_predict_iso
"""
    construct_predict_iso(df, cells, input_area="CA1", type=:has;
            dep_area=nothing, other_vars=[], other_ind_vars=[])

Returns a list of GLM formulae for predicting if a cycle contains
an isolated spike or the count of isolated spikes (but not neuron
specific). Just the count or has for the whole cycle.
"""
function construct_predict_iso(df, cells, input_area="CA1", type=:has;
        dep_area=nothing, other_vars=[], other_ind_vars=[])
    uArea = unique(cells.area)
    @assert length(uArea) == 2 "Only supports two area dataframes"
    dep_area = dep_area === nothing ? 
                setdiff(uArea, [input_area]) : dep_area
    ind_neurons = @subset(cells,:area .==input_area).unit
    filter!(n->string(n) ∈ names(df), ind_neurons)
    formulae = Vector{FormulaTerm}()
    ni = first(ind_neurons)
    independents = GLM.Term(Symbol(ni))
    for ni in ind_neurons[2:end]
        independents += GLM.Term(Symbol(ni)) 
    end
    if type == :has
        @info "type=$type"
        if dep_area == "CA1"
            formula = GLM.Term(:has_iso) ~ independents
        elseif dep_area == "PFC"
            formula = GLM.Term(:pfc_has_iso) ~ independents
        else
            throw(ErrorException("dep_area=$dep_area is unrecognized"))
        end
    elseif type == :count
        if dep_area == "CA1"
            formula = GLM.Term(:isolated_sum) ~ independents
        elseif dep_area == "PFC"
                formula = GLM.Term(:pfc_isolated_sum) ~ independents
        else
            throw(ErrorException("dep_area=$dep_area is unrecognized"))
        end
    else
        throw(ErrorException("$type is unrecognized"))
    end
    push!(formulae, formula)
    formulae
end

# -----------------
# Caching the dicts
# -----------------

export get_dx_dy
"""
    function get_dx_dy(df::DataFrame, relcyc::Int)

    Get the dataframes for the relative cycle `relcyc` and the cycle 0 data.

    # Arguments
    - `df::DataFrame`: The dataframe to get the data from.
    - `relcyc::Int`: The relative cycle to get the data from.

    # Returns
    - `dx::DataFrame`: The dataframe of the relative cycle `relcyc` data.
    - `dy::DataFrame`: The dataframe of the cycle 0 data.

    # Example
    ```julia
    df = DataFrame(x = 1:10, y = 1:10)
    df[:relcyc] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    dx, dy = get_dx_dy(df, 1)
    ```
    # Output
    ```
    (2×2 DataFrame
      Row │ x     y
          │ Int64 Int64
    ─────┼─────────────
        1 │     2     2
        2 │     3     3,
     2×2 DataFrame
      Row │ x     y
          │ Int64 Int64
    ─────┼─────────────
        1 │     1     1
        2 │     2     2)
    ```
"""
function get_dx_dy(df::DataFrame, relcyc::Int)
    @infiltrate
    dx = @subset(df, :relcycs .== relcyc)
    dy = @subset(df, :relcycs .== 0)
    _register_frames(dx, dy)
end

"""
         function _register_frames(dx::DataFrame, dy::DataFrame; register= [:cyc_batch, :cyc_match])

    Register the two dataframes `dx` and `dy` by the columns `register`.

    # Arguments
    - `dx::DataFrame`: The first dataframe to register.
    - `dy::DataFrame`: The second dataframe to register.
    - `register::Vector{Symbol}`: The columns to register by.

    # Returns
    - `dx::DataFrame`: The first dataframe with the columns `register` removed.
    - `dy::DataFrame`: The second dataframe with the columns `register` removed.

    # Example
    ```julia
    df = DataFrame(x = 1:10, y = 1:10)
    df[:cyc_batch] = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5]
    df[:cyc_match] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    df[:relcycs] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    df1 = @subset(df, :relcycs .== 1)
    df2 = @subset(df, :relcycs .== 0)
    dx, dy = _register_frames(df1, df2)
    ```
    # Output
    ```
    (2×2 DataFrame
      Row │ x     y
          │ Int64 Int64
    ─────┼─────────────
        1 │     2     2
        2 │     3     3,
     2×2 DataFrame
      Row │ x     y
          │ Int64 Int64
    ─────┼─────────────
        1 │     1     1
        2 │     2     2)
    ```
"""
function _register_frames(dx, dy; register= [:cyc_batch, :cyc_match])
    # Sort the dataframes by the register columns
    dx, dy = sort(dx, register), sort(dy, register)
    # Get the intersection of the register columns
    dxr, dyr = Matrix(dx[!,register]), Matrix(dy[!,register])
    # If the register columns are not the same, then get the intersection
    if dxr != dyr
        dxr, dyr = eachrow.((dxr, dyr))
        D = intersect(dxr, dyr)
        idx = [findfirst((x == d for x in dxr)) for d in D]
        idy = [findfirst((y == d for y in dyr)) for d in D]
        dx, dy = dx[idx,:], dy[idy, :]
    end
    # Assert that the register columns are the same
    @assert Matrix(dx[!,register]) == Matrix(dy[!,register])
    # Assert that the register columns are the same size
    @assert size(dx,1) == size(dy,1), "dx and dy are not the same size"
    dx, dy
end

"""
    get_futurepast_blocks(df)

Get the future and past blocks of the dataframe `df`.
"""
function get_futurepast_blocks(df)
    dxf = unstack(@subset(df, :relcycs .> 0), )
    dxp = @subset(df, :relcycs .<= 0)
    dy = @subset(df, :relcycs .== 0)
    dxp = groupby(dxp, [:cyc_batch, :cyc_match])
    dxf = groupby(dxf, [:cyc_batch, :cyc_match])
    dy = groupby(dy, [:cyc_batch, :cyc_match])
    kxp, kxf, ky = keys(dxp.keymap), keys(dxf.keymap), keys(dy.keymap)
    k = intersect(kxf,kxp,ky)
    dxf = sort(vcat([dxf[kk] for kk in k]...), [:cyc_batch, :cyc_match], )
    dxp = sort(vcat([dxp[kk] for kk in k]...), [:cyc_batch, :cyc_match], )
    dy  = sort(vcat([dy[kk] for kk in k]...), [:cyc_batch, :cyc_match], )
    Dict("pastcurr"=>(dxp, dy), "future"=>(dxf, dy))
end

function _handle_args(D::AbstractDict, area::String, 
    func::Function=x->adjr2(x, :devianceratio); kws...)
    _handle_args(Table.to_dataframe(D, func), area; kws...)
end
function _handle_args(D::DataFrame, area::String; 
        groupargs=nothing, combineargs=nothing, kws...) 
    if groupargs !== nothing
        @info "groupby" groupargs combineargs
        D = groupby(D, groupargs...) 
        D = combine(D, combineargs)
    end
    D, area
end

# -----------------
#   Plotting
# -----------------

@userplot GlmPlot
"""
    glmplot

takes either a dataframe of results or an abstract dict with a function to
instruct how to process the glm linear model objects
"""
@recipe function glmplot(plt::GlmPlot; groupargs=nothing, combineargs=nothing)
    D, area = _handle_args(plt.args...; groupargs, combineargs)
    data = sort(@subset(D, :indep .== area),:relcyc)
    @series begin
        seriestype := :hline
        linestyle := :dash
        c := :black
        label := ""
        ([0],)
    end
    @series begin
        seriestype := :vline
        linestyle := :dash
        c := :black
        label := ""
        ([0],)
    end
    alpha --> 0.5
    fillrange --> 0
    ylims --> (minimum(data.value), maximum(data.value))
    (data.relcyc, data.value)
end

# -----------------
# preparing cycles
# -----------------

export match_cycles!
"""
    match_cycles!(cycles::DataFrame, Rdf::DataFrame,
        occ::IndexedAdaptiveOcc; matches=3, iso_cycles = nothing)
find null cycles non-iso spike cycles matched on behavior per
iso spike cycle
# Arguments
- `cycles::DataFrame`: cycles dataframe
- `Rdf::DataFrame`: firing rate matrix R in the form of a dataframe
- `occ::IndexedAdaptiveOcc`: indexed adaptive occupancy object
- `matches::Int`: number of matches to find per cycle
- `iso_cycles::Vector{Int}`: vector of cycles containing isolated spikes
# Output
- `cycles::DataFrame`: cycles dataframe with a new column `matched` that
                       contains the matched cycles in terms of the behavior
                       as noted in the `occ` object
"""
function match_cycles!(cycles::DataFrame, Rdf::DataFrame, 
    occ::IndexedAdaptiveOcc; matches=3,
    iso_cycles = nothing, conditions=[:hasocc .== true])
    if iso_cycles === nothing
       throw(ArgumentError("iso_cycles must be provided"))
    end
    DIutils.filtreg.register(cycles, Rdf, transfer=["hasocc"], on="cycle")
    
    cycles.matched = Vector{Union{Vector{Int32}, Missing}}(missing,
        size(cycles,1))
    Threads.@threads for cyc in iso_cycles
        poss = [] 
        # Lookup cycles that match this isolated spike cycle's animal 
        # behavior
        ismissing(occ.datainds[cyc]) ? continue : nothing
        for gridmatch in occ.datainds[cyc]
            for cycmatch in occ.inds[gridmatch]
                push!(poss, cycmatch)
            end
        end
        # Lookup which cycles lack isolated spikes
        poss = poss[cycles[poss, :].isolated_sum .=== 0]
        samples_to_grab = min(length(poss), matches)
        if samples_to_grab > 0
            cycles[cyc,:].matched = 
                    Distributions.sample(poss, samples_to_grab, 
                    replace=false)
        end
    end
end

export df_FRpercycle_and_matched
"""
    df_FRpercycle_and_matched(cycles, Rdf_cycles, beh, val; iso_cycles,
        threading=true, indexers=[:time, :isolated_sum, :pfcisosum],
        cycrange=8)

obtain the fr per theta cycle of interest, relative cylces to it, and
cycles without isolated spikes matched on behavior

this grabs a span of cycles around each isolated spike cycle
and grabs a span around the null matched non-isolated cyles for each
iteration of the loop

# Arguments
- `cycles`: the cycles dataframe
- `Rdf_cycles`: dataframe with cyclewise groups of firing rate data
- `beh`: the behavior dataframe
- `firingvals`: values to mean (usually the column containing the firing
                rate)
- `iso_cycles`: the cycles with isolated spikes
- `threading`: whether to use threading
- `indexers`: the columns to use to index the cycles dataframe
- `cycrange`: the number of cycles to grab around the cycle of interest

# Return
- `df`: a vector of dataframes, each containing the firing rate data
- `cyc_error`: a dictionary of errors that occured during the loop
"""
function df_FRpercycle_and_matched(cycles, Rdf_cycles, beh, firingvals;
    iso_cycles, threading::Bool=true, indexers=[:time, :isolated_sum,
        :pfcisosum], cycrange=8)

    if threading
        nthreads = Threads.nthreads()
    else
        nthreads = 1
    end
    df = Vector{Vector{Union{Missing,DataFrame}}}(undef, nthreads)
    for thread in 1:nthreads
        df[thread] = Vector{Union{Missing,DataFrame}}()
    end
    # matched_cycle_holder = 
    #         Vector{Union{Int,Missing}}(missing, opt["matched"])
    cyc_error =  Dict() 
    Infiltrator.clear_disabled!()
    prog = Progress(length(iso_cycles); 
                     desc="grabing cycle batches into df")
    V = [firingvals, :i]
    E, M = Threads.Atomic{Int}(0), Threads.Atomic{Int}(0)

    # Grab a span of cycles around each isolated spike cycle
    # and grab a span around the null matched non-isolated cyles for each
    Threads.@threads for (i,cyc) in collect(enumerate(iso_cycles))
        # unit = parse(Int,replace(string(f.lhs), "_i"=>""))
        try
            tid = threading ? Threads.threadid() : 1
            cyc_batch = i
            # Push the isolated cycle and its preceding following cycles
            push!(df[tid], 
                grab_cycle_data(Rdf_cycles, cyc, V; indexers,
                                            cycrange=cycrange,
                                            cyc_central=cyc,
                                            cyc_batch, cyc_match=0))
            matched_cycs = @subset(cycles, :cycle .== cyc)
            matched_cycs = matched_cycs.matched[1]
            # Push MATCHED cycles
            if ismissing(matched_cycs) 
                Threads.atomic_add!(M, 1)
                @warn "Missing matches for $cyc, skipping"
                continue
            elseif isempty(matched_cycs)
                Threads.atomic_add!(M, 1)
                continue
            end
            for (j,mc) in enumerate(matched_cycs)
                push!(df[tid], 
                    grab_cycle_data(Rdf_cycles, mc, V; indexers, 
                                    cyc_central=mc,
                                    cycrange=cycrange,
                                    cyc_batch, cyc_match=j))
            end
            next!(prog)
        catch exception
            print("exception")
            cyc_error[cyc] = exception
            Threads.atomic_add!(E, 1)
        #     if mod(i, 100) == 0
        #         @info cyc_error
        #     end
            sleep(0.5)
            next!(prog)
        end
    end

    printstyled("Cycles without match ", M[]/length(iso_cycles), 
          "\nErrored cycles ", E[]/length(iso_cycles), color=:blink)

    # Because we are using threading, we need to combine the dataframes
    # we aggregated from each thread
    dfs = Vector{DataFrame}(undef, length(df))
    @info "Combining each thread's list"
    for i in eachindex(df)
        df[i]  = DIutils.arr.undef_to_missing(df[i])
        dfs[i] = vcat(skipmissing(df[i])...) # TODO this is slow
    end
    @info "Combining all threads' lists"
    df = vcat(dfs...)

    printstyled("Isocycles ", length(iso_cycles), 
        "\ntotal df cycles ", size(df,1),
        "\ntheoretical cycles ", length(iso_cycles) * ((2*cycrange)+1) * 4,
        color=:blink)
    println("")
    
    # Make sure all cycles have a match and that the number of matches
    # is greater than 1, otherwise we have a problem
    @assert :cyc_match ∈ propertynames(df) ||
        unique(df.cyc_match)>1 "FUCK"

    # Add column indicating if cycle has isolated spikes
    df.has_iso = df.isolated_sum .> 0
    # Obtain the columns who encode the firing rate of each neuron
    neuroncols = names(df)[tryparse.(Int, names(df)) .!== nothing]
    
    # Convert the firing rate values to the number of spikes per cycle
    df[:,neuroncols] .*= median(diff(beh.time)) 
    # Round the firing rate values to the nearest integer
    df[:,neuroncols] .= round.(df[:,neuroncols])
    # Convert the firing count values to integers
    df = DataFrames.transform(df, 
    neuroncols .=> n -> convert(Vector{Int64}, n), 
        renamecols=false)
    # Clean data frame, remove columns with all zeros
    col_all_zero = map(v->all(skipmissing(v.==0)), eachcol(df))
    df = df[!, Not(names(df)[col_all_zero])]

    return (;df, cyc_error)
end

export grab_cycle_data
"""
    grab_cycle_data(Rdf_cycles::GroupedDataFrame, cyc::Union{Int64,Int32}, 
        val::Symbol; indexers, cycrange::Int=8, kws...)

Takes a firing rate data frame grouped by cyles and grabs a span of
cycles around the cycle of interest.

Parameters
----------
- `Rdf_cycles::GroupedDataFrame`: Grouped data frame of firing rates
- `cyc::Union{Int64,Int32}`: Cycle of interest
- `val::Symbol`: Value to grab from the data frame
- `indexers::Vector{Symbol}`: Indexers to grab from the data frame
- `cycrange::Int=8`: Number of cycles to grab before and after the cycle
                     of interest 
- `kws...`: Keyword arguments to pass to `grab_cycle_data`

Returns
-------
- `DataFrame`: Data frame of the meaned firing rates per cycle spanning
                each cycle of interest
""" 
function grab_cycle_data(Rdf_cycles::GroupedDataFrame, 
        cyc::Union{Int64,Int32}, val::Symbol; indexers, 
        cycrange::Int=8, kws...)::DataFrame

     # Address cycles of interest
     🔑s = [(;cycle=cyc) 
            for cyc in UnitRange(cyc-cycrange, cyc+cycrange)
           ]

     # Which cycles to not compute the mean on in the unit unstacked
     # form of the data below
     do_not_mean = :area in propertynames(Rdf_cycles) ? 
                        Not([:time, :area]) : Not(:time)
    # Grab each cycle of activity
    U = [begin
            # TODO investigate nonunque
             u = unstack(Rdf_cycles[🔑], indexers, :unit, val,
                         combine=mean) 
             u = combine(u, do_not_mean .=> [mean], renamecols=false)
         end
        for 🔑 in 🔑s if 🔑 in keys(Rdf_cycles)]

    cycs = [🔑.cycle for 🔑 in 🔑s 
            if 🔑 in keys(Rdf_cycles)]
    relcycs = [🔑.cycle-cyc for 🔑 in 🔑s 
               if 🔑 in keys(Rdf_cycles)]

    # Added df to list
    df = DataFrames.hcat(DataFrame([cycs,relcycs],[:cycs,:relcycs]), 
                         vcat(U...; cols=:union))

    for (key, val) in kws
        df[!,key] .= val
    end
    
    return df
end

"""
    grab_cycle_data(Rdf_cycles::GroupedDataFrame, cyc::Union{Int64,Int32}, 
        vecofval::Vector{Symbol}; indexers, cycrange::Int=8, kws...)::DataFrame
    Grab a cycle of activity and its preceding and following cycles
    from a grouped data frame of cycles.

# Arguments
- `Rdf_cycles::GroupedDataFrame`: Grouped data frame of cycles
- `cyc::Union{Int64,Int32}`: Cycle of interest
- `vecofval::Vector{Symbol}`: Vector of values to grab per cycle
- `indexers::Symbol`: Indexer to grab per cycle
- `cycrange::Int=8`: Number of cycles to grab before and after
- `kws...`: Keyword arguments to add to the data frame

# Returns
- `DataFrame`: Data frame of the cycle of interest and its preceding and
"""
function grab_cycle_data(Rdf_cycles::GroupedDataFrame,
    cyc::Union{Int64,Int32}, vecofval::Vector{Symbol}; indexers, 
    cycrange::Int=8, kws...)::DataFrame
    dfs::Vector{DataFrame} = 
        [grab_cycle_data(Rdf_cycles, cyc, v; indexers, cycrange, kws...)
            for v in vecofval]
    out = dfs[1]
    on = ["cycs", "relcycs", "cyc_batch", "cyc_match"]
    mutualvars = names(out)[tryparse.(Int,names(out)) .=== 
                    nothing]
    mutualvars = union(mutualvars,on)
    for (v,i) in zip(vecofval[2:end], eachindex(dfs)[2:end])
        # df = leftjoin(df, dfs[i][!,Not(mutualvars)]; 
        #         on, makeunique=true, renamecols=""=>"i")
        @assert all([all(out[!,col] .== dfs[i][!,col])
            for col in on]) "Must match indices"
        out = hcat(out, dfs[i][!,Not(mutualvars)]; 
                    makeunique=true)
        renames = Dict(x => replace(x, "_1"=>"_$(string(v))")
                    for x in names(out) if occursin("_1",x))
        rename!(out, renames)
    end
    @assert "cyc_match" ∈ names(out)
    return out
end

"""
    ready_glm_vars(f::FormulaTerm, XX::DataFrame, y::DataFrame;
                   xtrans=identity, ytrans=identity, kws...)

Readies vars to be fed into the various glm functions 

# Arguments
- `f::FormulaTerm`: The formula to be used in the GLM
- `XX::DataFrame`: The data frame containing the predictors
- `y::DataFrame`: The data frame containing the response
- `xtrans::Function`: A function to be applied to the predictors
- `ytrans::Function`: A function to be applied to the response
- `expand_relcycs::Bool`: If true, will expand the relative cycles
- `zscoreX::Bool`: If true, will zscore the predictors
- `dummy_coding::Bool`: If true, will dummy code the predictors

# Returns
- `Tuple`: A tuple containing the predictors, response

# Example
```julia
using GLM, DataFrames, StatsModels
df = DataFrame(x=[1,2,3,4,5], y=[1,2,3,4,5])
f = @formula(y ~ x)
ready_glm_vars(f, df, df)
glm(f, df, df)
```
"""
function ready_glm_vars(f::FormulaTerm, XX::DataFrame, y::DataFrame;
xtrans=identity, ytrans=identity, expand_relcycs::Bool=false,
zscoreX::Bool=false, dummy_coding::Bool=false, kws...)::Tuple

    function get_mat(XX, f)
          cols = [string(ff) for ff in f.rhs]
          Matrix(XX[!,cols])
    end
    XX = if expand_relcycs && length(unique(XX.relcycs)) .> 1
        @debug "Expanding relative cycles"  
        unstack_relcycles(XX)
    else
        get_mat(XX, f)
    end
    y  = Vector(y[!,string(f.lhs)])
    @assert f.lhs ∉ f.rhs
    misses = (!).(ismissing.(y))
    XX, y = xtrans.(XX[misses,:]), ytrans.(y[misses])
    XX =  if dummy_coding  
        @debug "dummy coding"
        DIutils.statistic.dummycode(XX)
    elseif zscoreX
        @debug "zscoring"
        hcat([zscore(x) for x in eachcol(XX)]...)
    else
        XX
    end
    nansearch = isnan.(XX) .|| isnan.(y)
    goodcols = vec((!).(all(nansearch, dims=1)))
    goodrows = vec((!).(all(nansearch, dims=2)))
    XX, y = XX[goodrows, goodcols], y[goodrows]
    XX, y
end

"""
    expand_relcycles(XX::DataFrame, f::FormulaTerm)::DataFrame

Expands the relative cycles in the dataframe `XX`
"""
function unstack_relcycles(XX::DataFrame)::DataFrame
    res = []
    for cyc in sort(unique(XX.relcycs))
        push!(res, get_mat(@subset(XX, :relcycs .== cyc), f))
    end
    hcat(res...)
end

export expit
expit(x) = 1/(1+exp(-x))


"""
    measure_glm_dict!(D, y="y", ypred="ypred")

Adds metrics to the dictionary `D` and returns it.

# Inputs
- `D::Dict`: Dictionary to add metrics to
- `y::Symbol`: Name of the key for the true values
- `ypred::Symbol`: Name of the key for the predicted values

# Outputs
- `D::Dict`: Dictionary with metrics added
"""
function measure_glm_dict!(D::AbstractDict, y="y", ypred="ypred")
    y, ypred = D[y], D[ypred]
    M = Dict(
        "mae" => MLJ.mean_absolute_error(y, ypred),
        "cor" => cor(y, ypred),
        "r2" => MLJ.rsquared(y, ypred),
        "adjr2" => Metrics.adjusted_r2_score(y, ypred, length(y)),
    )
    merge!(D,M)
end


# -----------------------------------------
# The following structs are empty types
# used to dispatch on the various GLM methods
# -----------------------------------------

struct Custom_ElasticNetMLJ
end
struct SpecificGLM_MLJ
end
struct PYGLM_ElasticNet
end
struct GLMJL
end
struct GLM_MATLAB
end
struct GLM_R
end

GLMOUT = AbstractDict{String, Any}

export glm_
"""
    glm_(type, XX, y, Dist=Binomial(); kws...)

Wrapper for the various GLM methods

# Arguments
- `f::FormulaTerm`: The formula to be used in the GLM
- `XX::DataFrame`: The data frame containing the predictors
- `y::DataFrame`: The data frame containing the response
- `type`: The type of GLM to use. Can be one of the following:
    - `GLMJL`: Uses the GLM.jl package
    - `GLM_MATLAB`: Uses the GLM_MATLAB package
    - `GLM_R`: Uses the GLM_R package
    - `SpecificGLM_MLJ`: Uses the MLJ specific GLM
    - `Custom_ElasticNetMLJ`: Uses the MLJ specific ElasticNet
    - `PYGLM_ElasticNet`: Uses the Python specific ElasticNet

# Returns
- `Dict`: A dictionary containing the results of the GLM
"""
function glm_(f::FormulaTerm, XX::DataFrame, y::DataFrame, 
                 Dist=Binomial(); xtrans=identity, ytrans=identity,
                desc::AbstractDict=Dict(), kws...)
    @assert size(XX,1) == size(y,1) "Rows must match"
    XX, y = ready_glm_vars(f, XX, y; xtrans, ytrans, kws...)
    remove = [:zscoreX, :dummy_coding, :expand_relcycs]
    desc = merge(desc, Dict(k=>v for (k,v) in kws if k in remove))
    kws = DIutils.namedtup.pop(kws, remove)
    merge(glm_(XX, y, Dist; kws...), desc, Dict("formula"=>f))
end

export cvglm_
"""
    cvglm_(f::FormulaTerm, XX::DataFrame, y::DataFrame, 
                 Dist=Binomial(); xtrans=identity, ytrans=identity,
                desc::AbstractDict=Dict(), cv=5, shuffle=true, kws...)

 Cross-validated GLM that will call the wrapper to various GLM methods
`glm_` and return a dictionary containing the results of the 
 cross-validation
"""
function cvglm_(f::FormulaTerm, XX::DataFrame, y::DataFrame, 
                 Dist=Binomial(); xtrans=identity, ytrans=identity,
                desc::AbstractDict=Dict(), cv=5, shuffle=true, kws...)::GLMOUT
    @assert size(XX,1) == size(y,1) "Rows must match"
    XX, y = ready_glm_vars(f, XX, y; xtrans, ytrans, kws...)
    remove = [:zscoreX, :dummy_coding, :expand_relcycs]
    desc = merge(desc, Dict(k=>v for (k,v) in kws if k in remove))
    kws = DIutils.namedtup.pop(kws, remove)
    # Generate cross-validation sets using MLJ
    RES = OrdereredDict()
    cv = MLJ.CV(;nfolds=cv, shuffle)
    for (i,(train, test)) in enumerate(train_test_pairs(cv, 1:length(y)))
        XXtrain, ytrain = XX[train, :], y[train]
        XXtest,  ytest  = XX[test, :],  y[test]
        # Run the GLM
        res = glm_(XXtrain, ytrain, Dist; XXtest, ytest, kws...)
        merge!(res, desc, Dict("formula"=>f, "cv"=>i))
        push!(RES, (;cv=i)=>res)
    end
    RES
end

export shuffle_glm_
"""
    shuffle_glm_(f, XX, y, pos...; kwargs...)

Shuffle the cycle labels and run the GLM

# Arguments
- `f::FormulaTerm`: The formula to be used in the GLM
- `XX::DataFrame`: The data frame containing the predictors
- `y::DataFrame`: The data frame containing the response
- `pos...`: The positional arguments to pass to `glm_`
- `kwargs...`: The keyword arguments to pass to `glm_`

# Returns
- `Dict`: A dictionary containing the results of the GLM
"""
function shuffle_glm_(f, XX, y, pos...; shufflekws=Dict(), kwargs...)
    XX = shuffle_cyclelabels(XX; shufflekws...)
    glm_(f, XX, y, pos...; kwargs...)
end

export shuffle_cvglm_
"""
    shuffleglm_(f, XX, y, pos...; kwargs...)

Shuffle the cycle labels and run the GLM with cross-validation

# Arguments
- `f::FormulaTerm`: The formula to be used in the GLM
- `XX::DataFrame`: The data frame containing the predictors
- `y::DataFrame`: The data frame containing the response
- `pos...`: The positional arguments to pass to `glm_`
- `kwargs...`: The keyword arguments to pass to `glm_`

# Returns
- `Dict`: A dictionary containing the results of the GLM
"""
function shuffle_cvglm_(f, XX, y, pos...; shufflekws=Dict(), kwargs...)
    XX = shuffle_cyclelabels(XX; shufflekws...)
    cvglm_(f, XX, y, pos...; kwargs...)
end

"""
    glm_(F::Array{FormulaTerm}, pos...; kwargs...)

Iterate glm_ over an array of formula terms

# Arguments
- `F::Array{FormulaTerm}`: The array of formula terms
- `pos...`: The positional arguments to pass to `glm_`
- `kwargs...`: The keyword arguments to pass to `glm_`
- `handle_exception`: How to handle exceptions. Can be one of the following:
    - `:error`: Throw an error
    - `:warn`: Warn and continue
    - `:infiltrate`: Infiltrate and continue
    - `:exfiltrate`: Exfiltrate and continue
    - `:ignore`: Ignore and continue 

# Returns
- `Array{Dict}`: An array of dictionaries containing the results of the GLM
                for each formula term. The size of the array is the same
                as the size of `F`.
"""
function glm_(F::Array{FormulaTerm}, pos...; 
                handle_exception=:error, kwargs...)::GLMOUT
    out = []
    @showprogress for f in F
        try
        push!(out, glm_(f, pos...; kwargs...))
        catch exception
            if handle_exception == :warn
                @warn exception
            elseif handle_exception == :infiltrate
                @infiltrate
            elseif handle_exception == :exfiltrate
                @exfiltrate
            elseif handle_exception == :ignore
            else
                rethrow(exception)
            end
        end
    end
    reshape(out, size(F))
end

"""
    glm_(type, XX, y, Dist=Binomial(); kws...)

GLM function that dispatches on the type `type` to call the appropriate
GLM method.

# Inputs
- `XX::AbstractMatrix`: The predictors
- `y`: The response
- `Dist`: The distribution to use
- `type`: The type of GLM method to use

# Returns
- `Dict`: A dictionary containing the results
""" 
function glm_(XX::AbstractMatrix, y::AbstractVecOrMat, Dist=Binomial(); 
              ytest=nothing, type=:pyglm, kws...)

    ytest = ytest === nothing ? y : ytest

    # Select the type of GLM to dispatch to
    type = if type in [:mlj_spec, :spec]
         SpecificGLM_MLJ()
    elseif type in [:mlj_custom, :custom]
         Custom_ElasticNetMLJ()
    elseif type in [:pyglm, :pyglmnet]
         PYGLM_ElasticNet()
    elseif type == :matlab
         GLM_MATLAB()
    elseif type in [:r,:R]
         GLM_R()
    else
        throw(ArgumentError("type=$type not recognized"))
    end

    # Dispatch to the appropriate GLM method
    kv  = glm_(type, XX, ytest, Dist; kws...)

    # Merge the results with top level results and the values
    # to be predicted
    out = merge(Dict("type"=>string(type), "y"=>ytest,
                      "yr"=>denserank(ytest)./length(ytest), 
                      "ypredr"=>denserank(kv["ypred"])./length(ytest)
                ), kv)
    out["ypred"], out["y"] = Float64.(out["ypred"]), Float64.(out["y"])
    measure_glm_dict!(out)
    measure_glm_dict!(out, "yr", "ypredr")
end

"""
    glm_(::SpecificGLM_MLJ, XX, y, Dist=Binomial(); kws...)

Wrapper for the specific GLM method

# Inputs
- `XX::AbstractMatrix`: Design matrix
- `y::AbstractVector`: Response vector
- `Dist::Distribution`: Distribution to use
"""
function glm_(::PYGLM_ElasticNet, XX::AbstractMatrix,
y::Union{AbstractVector, AbstractMatrix}, Dist=Binomial();
testXX=nothing, ytest=nothing, kws...)
    ytest  = ytest === nothing ? y : ytest
    testXX = testXX === nothing ? XX : testXX
    kws = (;tol=1e-3, score_metric="pseudo_R2", alpha=0.5, 
            max_iter=100, cv=3, kws...)
    gl_glm = pyglmnet.GLMCV(;distr=diststr(Dist), kws...)
    gl_glm.fit(XX, y)
    ypred = gl_glm.predict(testXX)
    score = gl_glm.score(testXX, ytest)
    out = Dict(string.(collect(keys(kws))) .=> collect(values(kws)))
    merge(out, Dict("y"=>y, "ypred"=>ypred), 
        Dict("lambda"=>gl_glm.reg_lambda_opt_, "lambdas"=>gl_glm.reg_lambda,
        kws.score_metric=>score, "score"=>score)
    )
end
glm_pyglm(XX, y, Dist, kws...) = glm_(PYGLM_ElasticNet(), XX, y, Dist;
                                 kws...)

"""
    glm_(::Custom_ElasticNetMLJ, XX, y, Dist=Binomial(); kws...)

Wrapper for the custom GLM method

# Inputs
- `XX::AbstractMatrix`: Design matrix
- `y::AbstractVector`: Response vector
- `Dist::Distribution`: Distribution to use

# Outputs
Dictionary of results
"""
function glm_(::Custom_ElasticNetMLJ, XX::AbstractMatrix, y, Dist=Binomial();
    testXX=nothing, kws...)
    # testXX = testXX === nothing ? XX : testXX
    # @debug "Custom $(typeof(Dist))"
    # yin, XXin = if Dist isa Binomial || Dist isa Poisson
    #     replace(y, 0=>1e-16), replace(XX, 0=>1e-16)
    # else
    #     y
    # end
    # m = ypred = nothing
    # R = ElasticNetCVRegressor(n_jobs=Threads.nthreads(),
    #     cv=min(10,size(XX,1)), normalize=true)
    # lif(x) = GLM.linkfun.(canonicallink(Dist), x)
    # ilif(x) = GLM.linkinv.(canonicallink(Dist), x)
    # try
    #     XXin, yin = XX, lif(yin)
    #     XXin, y = MLJ.table(XXin), y
    #     m = MLJ.machine(R, XXin, yin)
    #     MLJ.fit!(m)
    #     ypred = ilif(MLJ.predict(m, XXin))
    # catch exception
    #     println("Exception = $exception")
    # end
    # Dict("m"=>m, "ypred"=>ypred, "coef"=> m.fitresult[1].coef_)
end
glm_custom(XX::AbstractMatrix, y::AbstractVector, Dist) = 
        glm_(Custom_ElasticNetMLJ(), XX::AbstractMatrix, y, Dist=Binomial())

"""
glm_(::SpecificGLM_MLJ, XX::AbstractMatrix, y, Dist)

mlj version of glmnet

# Inputs
- `XX::AbstractMatrix`: Design matrix
- `y::AbstractVector`: Response vector
- `Dist`: Distribution
"""
function glm_(::SpecificGLM_MLJ, XX::AbstractMatrix, y, Dist; 
             testXX=nothing, kws...)
    # testXX = testXX === nothing ? XX : testXX
    # if Dist isa Binomial
    #     @debug "Spec binoomial"
    #     XX, y = MLJ.table(XX), y
    #     R = LogisticCVClassifier(n_jobs=Threads.nthreads(),
    #         penalty="elastic_net",
    #         cv=min(5,size(XX[1],1)))
    #     m = MLJ.machine(R, XX, η)
    #     MLJ.fit!(m)
    #     ypred = MLJ.predict(m, XX)
    #     coef = m.fitresult[1]
    # elseif Dist isa Poisson
    #     @debug "Spec poisson"
    #     sklearn = pyimport("sklearn")
    #     lm = sklearn.linear_model
    #     R = lm.PoissonRegressor()
    #     R.fit(XX, y)
    #     ypred = R.predict(XX)
    #     m = nothing
    #     coef = R.coef_
    # end
    # Dict("m"=>m, "ypred"=>ypred, "coef"=>coef)
end
glm_specific(XX::AbstractMatrix, y, Dist) = 
        glm_(SpecificGLM_MLJ(), XX::AbstractMatrix, y, Dist=Binomial())

export glm_
"""
    glm_(::GLMJL, XX, y, Dist=Binomial(); testXX, kws...)

Wrapper for the specific GLM method from GLM.jl
"""
function glm_(::GLMJL, XX, y; testXX=nothing, kws...)
    testXX = testXX === nothing ? XX : testXX
    @debug "GLM.jl"
    y =  expit.(Int.(y[misses]) .> 0)
    m =  GLM.glm(XX, y, Dist)
    ypred = GLM.predict(XX)
    Dict("m"=>m, "ypred"=>ypred, "coef"=>GLM.coef(m))
end
glm_glmjl(XX::AbstractMatrix, y) = glm_(GLMJL(), XX, y)

# ------------------------
# METHOD : MATLAB
# `fitglm` or better `lassoglm`
# ------------------------

export glm_matlab
"""
    glm_matlab(f::FormulaTerm, XX::DataFrame, y::DataFrame,
                Dist=Binomial(); testXX=nothing, quick_and_dirty=false,
                pre=nothing, kws...)

Wrapper for the specific GLM method from matlab (fitglm and lassoglm)

# Inputs
- `XX::AbstractMatrix`: Design matrix
- `y::AbstractVector`: Response vector
- `Dist::Distribution`: Distribution to use
- `link::Link`: Link function to use
- `testXX::AbstractMatrix`: Test design matrix
- `quick_and_dirty::Bool`: Use `glmfit` instead of `lassoglm`
- `pre::Dict`: Precomputed results; if passed, will skip the fitting
"""
function glm_(::GLM_MATLAB, XX, y, Dist=Binomial(), 
    link=canonicallink(Dist); testXX=nothing,
    quick_and_dirty=false, pre=nothing, kws...)

    testXX = testXX === nothing ? XX : testXX

    dist = diststr(Dist)
    link = linkstr(link)

    @debug "matlab" link dist

    xnz, y = DIutils.arr.atleast2d(XX), DIutils.arr.atleast2d(y)
    @debug quick_and_dirty
    if pre === nothing
        if quick_and_dirty
           @time mat"[$B, ~, $stats] = glmfit(double($xnz), double($y), $dist);"
        else
           @time mat"[$B, $stats] = lassoglm(double($xnz), double($y), $dist, 'alpha', 0.5, 'CV', 3, 'MCReps', 5);"
        end
       y   = vec(Float64.(y))
       pre = Dict("B"=>B, "stats"=>stats, "type"=>:matlab, "y"=>y)
    end
    B, stats = pre["B"], pre["stats"]
    if quick_and_dirty
        c = B
    else
       idxLambdaMinDeviance = stats["IndexMinDeviance"]
       intercept = stats["Intercept"]
       B0   = intercept[Int(idxLambdaMinDeviance)];  
       c = [B0; B[ :,Int(idxLambdaMinDeviance) ]]
    end
    ypred = mat"glmval($c, double($(xnz)), $link, $stats);"
    y, ypred = vec(Float64.(y)), vec(Float64.(ypred))
    merge!(pre, Dict( "ypred"=>ypred, "y"=>y, 
         "coef"=>c[2:end], "intercept"=>c[1]))
    merge!(pre, pop!(pre, "stats"))
end
glm_matlab(XX::AbstractMatrix, y, Dist, kws...) =
        glm_(GLM_MATLAB(), XX, y, Dist; kws...)


# ------------------------
# METHOD : R
# ------------------------

using RCall
"""
    glm_(::GLM_R, XX::AbstractMatrix, y, Dist)

The R programming language version of the GLM    
"""
function glm_(::GLM_R, XX::AbstractMatrix, y, Dist; testXX=nothing, 
                  kws...)
    testXX = testXX === nothing ? XX : testXX
    # Import GLM net R libraries
    R"library(glmnet)"
    # Run GLM on XX and y with distribution Dist
    Dist = diststr(Dist)
    R"glmnet.fit <- glmnet($testXX, $y, family = $Dist, alpha = 0.5)"
    # Get the coefficients
    coef = R"glmnet.fit$beta"
    # Get the predictions
    ypred = rcopy(R"predict(glmnet.fit, $XX)")
    # Get the intercept
    intercept = rcopy(R"glmnet.fit$a0")
    # Get the deviance
    deviance = rcopy(R"glmnet.fit$dev.ratio")
    # Get the lambda
    lambda = rcopy(R"glmnet.fit$lambda")
    # Get the null deviance
    nulldeviance = rcopy(R"glmnet.fit$nulldev")
    # Get the number of iterations
    niter = rcopy(R"glmnet.fit$niter")
    # Return the results
    Dict("ypred"=>ypred, "intercept"=>intercept,
    "deviance"=>deviance, "lambda"=>lambda, "nulldeviance"=>nulldeviance,
    "niter"=>niter)
end


export run_glm!
"""
Run the GLM on several batches of data

- `df` : the dataframe, with the data
- `df_cache` : the dataframe, with the data
- `cells` : dataframe describing the possible cells
- `formula_method` : the function to create the formulae
- `glmtool` : the tool to use
- `Dist` : the distribution to use
- `unitwise` : if true, we do the GLM on the unitwise data or instead
               a different column specified by the formula_method
- `unitrep` : unit replace ie, when a unit is address with a string,
              a replacement argment to address something else (e.g. "_i")
             to pull an isolated spike column for that cell
- `xtrans` : the function to transform the X values
- `ytrans` : the function to transform the Y values
- `cacheXY` : a dictionary to cache the X and Y values
- `modelz` : a dictionary to cache the model values

# Returns
- `modelz` : a dictionary to cache the model values
- `cacheXY` : a dictionary to cache the X and Y values

The outputs are also written to by reference, so no need to return them.
"""
function run_glm!(full_df::DataFrame, df_batches::AbstractDict,
    cells::DataFrame, formula_method::Function, glmtool::Symbol; 
    area_to_area::Bool=false,
    Dist=Binomial(), unitwise::Bool, unitrep=[], xtrans::Function=identity,
    ytrans::Function=identity, cacheXY=ThreadSafeDict(),
    dryrun::Bool=false, modelz=ThreadSafeDict())
    # Now let's take the formula and apply them
    formulae = OrderedDict() 
    formulae["ca1pfc"] = formula_method(full_df, cells, "CA1");
    formulae["pfcca1"] = formula_method(full_df, cells, "PFC");
    if area_to_area
        formulae["ca1ca1"] = formula_method(full_df, cells, "CA1", 
                                     dep_area="CA1");
        formulae["pfcpfc"] = formula_method(full_df, cells, "PFC"; 
                                     dep_area="PFC");
    end
    @assert !isempty(first(values(formulae)))
    glmsets = []
    # CREATE THE ITERATOR
    for indep in keys(formulae), cachekey in keys(df_batches), 
       f in formulae[indep]
        push!(glmsets, (indep, cachekey, f))
    end
    # DRY RUN?
    # if dryrun, then we only do the first two
    if dryrun
        glmsets = glmsets[1:2]
    end
    prog = Progress(length(glmsets); desc="GLM spike counts")
    for (indep, x_key, f) in glmsets                                                

        # Get 🔑 the output key
        if unitwise
            unitstr = !isempty(unitrep) ?
                replace(string(f.lhs),unitrep...) : string(f.lhs)
            unit=parse(Int, unitstr)
            out_key = (;unit, indep, x_key...)
        else
            out_key = (;indep, x_key...)
        end

        # Run GLM model 🤖
        try
            XX, y = df_batches[x_key]
            modelz[(;indep, x_key, unit)] = glm_(f, XX, y, Dist;
                                   type=glmtool, xtrans=xtrans, ytrans=ytrans)
            cacheXY[out_key] = (;XX, y)
        catch exception
            modelz[(;indep, x_key, unit)] = exception
            print("Error on $out_key, exception = $exception")
            sleep(0.1)
        end
        next!(prog)
     end
    return modelz, cacheXY
end
"""
Runs the GLM on the full dataframe
"""
function run_glm!(full_df::DataFrame,
cells::DataFrame,formula_method::Function, glmtool::Symbol; kws...)
    df_batches = OrderedDict(
        (;batch=:full) => (full_df, full_df)
    )   
    run_glm!(full_df, df_batches, cells, formula_method, glmtool; kws...)
end
"""
Runs the GLM on a grouped dataframe
"""
function run_glm!(grouped::GroupedDataFrame, cells::DataFrame,
formula_method::Function, glmtool::Symbol; kws...)
    df_batches = OrderedDict(
        begin
            k = NamedTuple(grouped.cols, keytup)
            k => (grouped[keytup], grouped[keytup])
        end 
    for keytup in keys(grouped)
    )   
    run_glm!(DataFrame(grouped), df_batches, cells, formula_method,
                glmtool; kws...)
end


export shuffle_cyclelabels
"""
    shuffle_cyclelabels(df, perm=[:cycs,:relcycs])

Shuffle the cycle labels in the dataframe `df` and return the shuffled
dataframe.

# Arguments
 - `df`: a dataframe with columns `perm`
 - `perm`: a vector of symbols, the columns to shuffle

# Returns
 - `dfshuf`: a dataframe with the same columns as `df`, but with the
columns `perm` shuffled.
"""
function shuffle_cyclelabels(df, perm=[:cycs,:relcycs])
    # if "orig"*string(perm[1]) ∉ names(df)
    #         df[!,"orig".*string(perm)] .= df[!, perm]
    # end
    inds = collect(1:size(df,1))
    randperm!(inds)
    dfshuf = DataFrames.transform(df, perm .=> 
            x->x[inds], 
            Not(perm), renamecols=false)
    dfshuf = sort(dfshuf, perm)
    @assert dfshuf != df
    @assert dfshuf[!,perm] != df[!,perm]
    return dfshuf
end

# FOR INTERACTING WITH GLM MODEL DICTS
"""
    Table.to_dataframe(D::AbstractDict)::DataFrame

Convert a dictionary to a dataframe, where the keys are the column names
"""
function Table.to_dataframe(D::AbstractDict, func::Function)::DataFrame
    kt = keytype(D)
    if all(isa.(collect(values(D)), AbstractDict))
        D = Dict{kt, Any}(k=>Table.to_dataframe(v, func) for (k,v) in D)
        Table.to_dataframe(D)
    else
        D = Dict{kt, Any}(k=>func(k, v) for (k,v) in D)
        Table.to_dataframe(D)
    end
end


export clean_keys
"""
Clean the keys of a dictionary, removing the blacklist
"""
function clean_keys(D::AbstractDict, blacklist)
    Dict(k=>clean_keys(v, blacklist) for (k,v) in D if k ∉ blacklist) 
end
clean_keys(D::Any, blacklist) = D

export grabfield
grabfield(x::String) = (k,v)->if v isa Exception || v === nothing
    NaN
else
    try v[x]
    catch NaN end
end
grabfield(x::Vector{String}) = (k,v)->if v isa Exception || v === nothing
    Dict(k=>NaN for k in x)
else
    try Dict(k=>v[k] for k in x)
    catch; Dict(k=>NaN for k in x); end
end

export ucellcols
"""
grab columsn that correspond to cell props
"""
ucellcols(df) = names(df)[tryparse.(Int, names(df)) .!== nothing]
cellcols(df) = ucellcols(df)
export cellcols

export uicellcols
"""
grab columsn that correspond to cell isolated props
"""
function uicellcols(df)
    nms = [nm for nm in names(df) if endswith(nm,"_i")]
    names_cleaned = replace.(nms, ["_i"=>""])
    nms[tryparse.(Int, names_cleaned) .!== nothing]
end
isolatedcellcols(df) = uicellcols(df)
export isolatedcellcols


# STRING METHODS
nameit(x) = lowercase(begin
    first(split(first(split(string(x),"{")), "(" ))
end)

"""
    linkstr(link::Distribution)

Returns the name of the link function
"""
linkstr(Dist) = 
    if Dist isa Distribution
        replace(nameit(canonicallink(Dist)),"link"=>"")
    elseif Dist isa String
            Dist
    else
        replace(nameit(Dist),"link"=>"")
    end

export diststr
"""
    diststr(Dist::Distribution)

Returns the name of the distribution
"""
diststr(Dist::Union{String,Distribution}) = Dist isa String ? 
                Dist : nameit(Dist)


end
