__precompile__(false)
"""
    isolated

This module contains code used to analyze isolated spikes. It is
designed to be used in conjunction with the `GoalFetchAnalysis` module
and the `GoalFetchAnalysis.Munge` module.
"""
module isolated
    using MATLAB: matlab_cmd
    using ArgParse: test_range
    import DIutils
    using DIutils.binning
    import Logging

    using JLD2, ArgParse, DrWatson, DIutils.dict, RecipesBase, DataFramesMeta,
          Statistics, NaNStatistics, Infiltrator, DataFrames, 
            Missings, Plots, StatsBase, ThreadSafeDicts
    using DataStructures: OrderedDict
    import Distributions
    import DIutils: Table
    using ProgressMeter
    
    using GLMNet, MultivariateStats, MLJ, ScikitLearn, Metrics, GLM
    using MATLAB, PyCall
    using Random
    using MLJScikitLearnInterface: ElasticNetCVRegressor
    function __init__()
        pyglmnet = pyimport("pyglmnet")
    end

    init_mlj = false

    function __init__mlj()
        @eval isolated init_mlj = true
        if !init_mlj
            @eval isolated using MLJScikitLearnInterface: 
                 ElasticNetCVRegressor, LogisticCVClassifier
        end
    end

    export path_iso
    function path_iso(animal::String, day::Int, tet=:ca1ref)::String
        datadir("isolated","iso_animal=$(animal)_day=$(day)_tet=$(tet).jld2")
    end
    function path_iso(opt::AbstractDict)::String
        path_iso(opt["animal"], opt["day"], opt["tet"])
    end
    function path_iso(pos...; append::String)::String
        f = path_iso(pos...)
        replace(f, ".jld2" => "$(append).jld2")
    end

    export load_iso
    """
        load_iso(pos...)

    Returns a dictionary of the isolated variables of interest
    
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
    - `formulae`: a list of GLM formulae
    """
    function construct_predict_isospikecount(df, cells, input_area="CA1";
            other_vars=[], other_ind_vars=[])
        uArea = unique(cells.area)
        @assert length(uArea) == 2 "Only supports two area dataframes"
        dependent_area = setdiff(uArea, [input_area])

        dep_neurons = @subset(cells,:area .==dependent_area).unit
        ind_neurons = @subset(cells,:area .==input_area).unit
        dep_neurons = string.(dep_neurons) .* "_i"
        ind_neurons = string.(ind_neurons) 
        filter!(n->n ∈ names(df), dep_neurons)
        filter!(n->n ∈ names(df), ind_neurons)
        
        formulae = []
        for nd in dep_neurons
            ni = first(ind_neurons)
            independents = GLM.Term(Symbol(ni))
            for ni in ind_neurons[2:end]
                independents += GLM.Term(Symbol(ni)) 
            end
            formula = GLM.Term(Symbol(nd)) ~ independents
            push!(formulae, formula)
        end
        formulae
    end

    export construct_predict_spikecount
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
        
        formulae = []
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
    function construct_predict_iso(df, cells, input_area="CA1", type=:has;
            other_vars=[], other_ind_vars=[])
        uArea = unique(cells.area)
        @assert length(uArea) == 2 "Only supports two area dataframes"
        ind_neurons = @subset(cells,:area .==input_area).unit
        filter!(n->string(n) ∈ names(df), ind_neurons)
        formulae = []
        ni = first(ind_neurons)
        independents = GLM.Term(Symbol(ni))
        for ni in ind_neurons[2:end]
            independents += GLM.Term(Symbol(ni)) 
        end
        if type == :has
            @info "type=$type"
            formula = GLM.Term(:has_iso) ~ independents
        elseif type == :count
            formula = GLM.Term(:isolated_sum) ~ independents
        else
            throw(ErrorException("$type is unrecognized"))
        end
        push!(formulae, formula)
        formulae
    end

    # -----------------
    # Caching the dicts
    # -----------------

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
        dx, dy = sort(dx, register), sort(dy, register)
        dxr, dyr = Matrix(dx[!,register]), Matrix(dy[!,register])
        if dxr != dyr
            dxr, dyr = eachrow.((dxr, dyr))
            D = intersect(dxr, dyr)
            idx = [findfirst((x == d for x in dxr)) for d in D]
            idy = [findfirst((y == d for y in dyr)) for d in D]
            dx, dy = dx[idx,:], dy[idy, :]
        end
        @assert Matrix(dx[!,register]) == Matrix(dy[!,register])
        dx, dy
    end

    """
        get_futurepast_blocks(df)
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
    """match_cycles!
    
    find null cycles non-iso spike cycles matched on behavior per
    iso spike cycle
    """
    function match_cycles!(cycles::DataFrame, Rdf::DataFrame, 
        occ::IndexedAdaptiveOcc; matches=3,
        iso_cycles = nothing)

        cycles.hasocc = (!).(ismissing.(occ.datainds))

        if iso_cycles === nothing
           unique(@subset(Rdf, :isolated_sum .> 0, 
                          :hasocc .== true).cycle) 
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

    # vcatnonmiss(df::DataFrame) = DataFrames.vcat(df[(!).(ismissing.(df))]...)
    
    export df_FRpercycle_and_matched
    """
    df_FRpercycle_and_matches
    
    obtain the fr per theta cycle of interest, relative cylces to it, and
    cycles without isolated spikes matched on behavior

    changes
    """
    function df_FRpercycle_and_matched(cycles, Rdf_cycles, beh, val;
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
        V = [val, :i]
        E, M = Threads.Atomic{Int}(0), Threads.Atomic{Int}(0)
        #[(length(iso_cycles)-100):end]
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
                sleep(0.05)
                next!(prog)
            end
        end
        printstyled("Cycles without match ", M[]/length(iso_cycles), 
              "\nErrored cycles ", E[]/length(iso_cycles), color=:blink)
        dfs = Vector{DataFrame}(undef, length(df))
        for i in eachindex(df)
            dfs[i] = vcat(skipmissing(df[i])...)
        end
        df = vcat(dfs...)
        printstyled("Isocycles ", length(iso_cycles), 
            "\ntotal df cycles ", size(df,1),
            "\ntheoretical cycles ", length(iso_cycles) * ((2*cycrange)+1) * 4,
            color=:blink)
        @assert :cyc_match ∈ propertynames(df) ||
            unique(df.cyc_match)>1 "FUCK"
        df.has_iso = df.isolated_sum .> 0
        # Spike count
        neuroncols = names(df)[tryparse.(Int, names(df)) .!== nothing]
        # TODO not INT because it's gaussian smoothed
        df[:,neuroncols] .*= median(diff(beh.time)) 
        df[:,neuroncols] .= round.(df[:,neuroncols])
        df = DataFrames.transform(df, 
        neuroncols .=> n -> convert(Vector{Int64}, n), 
            renamecols=false)
        # Clean data frame
        col_all_zero = map(v->all(skipmissing(v.==0)), eachcol(df))
        df = df[!, Not(names(df)[col_all_zero])]
        (;df, cyc_error)
    end

    export grab_cycle_data
    function grab_cycle_data(Rdf_cycles::GroupedDataFrame, 
            cyc::Union{Int64,Int32}, val::Symbol; indexers, 
            cycrange::Int=8, kws...)::DataFrame
         selector = :area in propertynames(Rdf_cycles) ? Not([:time, :area]) : 
                                                         Not(:time)
         # Address cycles of interest
         🔑s = [(;cycle=cyc) 
                for cyc in UnitRange(cyc-cycrange, cyc+cycrange)
               ]

        # Grab each cycle of activity
        U = [begin
                # TODO investigate nonunque
                 u = unstack(Rdf_cycles[🔑], indexers, :unit, val,
                             combine=last) 
                 u = combine(u, selector .=> [mean], renamecols=false)
             end
            for 🔑 in 🔑s if 🔑 in keys(Rdf_cycles)]
         # @info combine(groupby(Rdf_cycles[🔑],:unit),
         #               :time=>x->length(x)==length(unique(x)))

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
        
        df
    end
    
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
        out
    end

    """
        ready_glm_vars(f::FormulaTerm, XX::DataFrame, y::DataFrame;
                       xtrans=identity, ytrans=identity, kws...)

    Readies vars to be fed into the various glm functions 
    """
    function ready_glm_vars(f::FormulaTerm, XX::DataFrame, y::DataFrame;
                    predictkey=0, xtrans=identity, ytrans=identity, 
                    zscoreX=false, dummy_coding=false, kws...)::Tuple

        function get_mat(XX, f)
              cols = [string(ff) for ff in f.rhs]
              Matrix(XX[!,cols])
        end
        XX = if length(unique(XX.relcycs)) .> 1
            res = []
            for cyc in sort(unique(XX.relcycs))
                push!(res, get_mat(@subset(XX, :relcycs .== cyc), f))
            end
            hcat(XX...)
        else
            get_mat(XX, f)
        end
        y  = Vector(y[!,string(f.lhs)])
        @assert f.lhs ∉ f.rhs
        misses = (!).(ismissing.(y))
        XX, y = xtrans.(XX[misses,:]), ytrans.(Int.(y[misses]))
        XX =  if dummy_coding  
                @info "dummy coding"
                DIutils.statistic.dummycode(XX)
            elseif zscoreX
                @info "zscoring"
                hcat([zscore(x) for x in eachcol(XX)]...)
            else
                XX
            end
        XX, y
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


    # -------------------------
    # METHOD : MLJ package GLM  
    # Which taps python's sklearn
    # Also added another method
    # -------------------------

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

    export glm_
    """
        glm_(type, XX, y, Dist=Binomial(); kws...)
    
    Wrapper for the various GLM methods
    """
    function glm_(f::FormulaTerm, XX::DataFrame, y::DataFrame, 
                     Dist=Binomial(); xtrans=identity, ytrans=identity,
                    kws...)
        @assert size(XX,1) == size(y,1) "Rows must match"
        XX, y = ready_glm_vars(f, XX, y; xtrans, ytrans, kws...)
        kws = DIutils.namedtup.pop(kws, [:zscoreX, :dummy_coding])
        glm_(XX, y, Dist; kws...)
    end

    function glm_(XX::AbstractMatrix, y, Dist=Binomial(); type=:specific, kws...)
         __init__mlj()
         type = if type == :mlj_spec
             SpecificGLM_MLJ()
         elseif type == :mlj_custom
             Custom_ElasticNetMLJ()
         elseif type == :pyglm
             PYGLM_ElasticNet()
         elseif type == :matlab
             GLM_MATLAB()
         else
            @error "type=$type not recognized"
         end
         kv  = glm_(type, XX, y, Dist; kws...)
        out = merge(Dict("type"=>string(type), "y"=>y,
                          "yr"=>denserank(y)./length(y), 
                          "ypredr"=>denserank(kv["ypred"])./length(y)
                    ), kv)
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
    function glm_(::PYGLM_ElasticNet, XX::AbstractMatrix, y, Dist=Binomial(); kws...)
        kws = (;tol=1e-3, score_metric="pseudo_R2", alpha=0.5, 
                max_iter=100, cv=3, kws...)
        @infiltrate
        gl_glm = pyglmnet.GLMCV(;distr=diststr(Dist), kws...)
        gl_glm.fit(XX, y)
        ypred = gl_glm.predict(XX)
        out = Dict(string.(collect(keys(kws))) .=> collect(values(kws)))
        merge(out, Dict("y"=>y, "ypred"=>ypred), 
        Dict("lambda"=>gl_glm.reg_lambda_opt_, "lambdas"=>gl_glm.reg_lambda))
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
    function glm_(::Custom_ElasticNetMLJ, XX::AbstractMatrix, y, Dist=Binomial())
        @info "Custom $(typeof(Dist))"
        yin, XXin = if Dist isa Binomial || Dist isa Poisson
            replace(y, 0=>1e-16), replace(XX, 0=>1e-16)
        else
            y
        end
        m = ypred = nothing
        R = ElasticNetCVRegressor(n_jobs=Threads.nthreads(),
            cv=min(10,size(XX,1)), normalize=true)
        lif(x) = GLM.linkfun.(canonicallink(Dist), x)
        ilif(x) = GLM.linkinv.(canonicallink(Dist), x)
        try
            XXin, yin = XX, lif(yin)
            XXin, y = MLJ.table(XXin), y
            m = MLJ.machine(R, XXin, yin)
            MLJ.fit!(m)
            ypred = ilif(MLJ.predict(m, XXin))
        catch exception
            println("Exception = $exception")
            @infiltrate
        end
        Dict("m"=>m, "ypred"=>ypred, "coef"=> m.fitresult[1].coef_)
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
    function glm_(::SpecificGLM_MLJ, XX::AbstractMatrix, y, Dist)
        if Dist isa Binomial
            @info "Spec binoomial"
            XX, y = MLJ.table(XX), y
            R = LogisticCVClassifier(n_jobs=Threads.nthreads(),
                penalty="elastic_net",
                cv=min(5,size(XX[1],1)))
            m = MLJ.machine(R, XX, η)
            MLJ.fit!(m)
            ypred = MLJ.predict(m, XX)
            coef = m.fitresult[1]
        elseif Dist isa Poisson
            @info "Spec poisson"
            sklearn = pyimport("sklearn")
            lm = sklearn.linear_model
            R = lm.PoissonRegressor()
            R.fit(XX, y)
            ypred = R.predict(XX)
            m = nothing
            coef = R.coef_
        end
        Dict("m"=>m, "ypred"=>ypred, "coef"=>coef)
    end
    glm_specific(XX::AbstractMatrix, y, Dist) = 
            glm_(SpecificGLM_MLJ(), XX::AbstractMatrix, y, Dist=Binomial())
    
    export glm_
    function glm_(::GLMJL, XX, y)
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
                 Dist=Binomial(); xtrans=identity, ytrans=identity,
                kws...)
    """
    function glm_(::GLM_MATLAB, XX, y, Dist=Binomial(), link=nothing; 
                        quick_and_dirty=false, pre=nothing, kws...)

        dist = diststr(Dist)
        link = linkstr(link)

        @info "matlab" link dist

        xnz, y = DIutils.arr.atleast2d(xnz), DIutils.arr.atleast2d(y)
        @info quick_and_dirty
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
        # ypred = mat"glmval($c, double($(xnz)), 'log', $stats)"
        #mat"[$ypred, $dlo, $dhi] = glmval($c, double($(XX)), 'log', $stats)"
        ypred = mat"glmval($c, double($(xnz)), $link, $stats);"
        y, ypred = vec(Float64.(y)), vec(Float64.(ypred))
        merge!(pre, Dict( "ypred"=>ypred, "y"=>y, 
             "coef"=>c[2:end], "intercept"=>c[1]))
    end
    glm_matlab(XX::AbstractMatrix, y, Dist, kws...) =
            glm_(GLM_MATLAB(), XX, y, Dist; kws...)


    export run_glm!
    """
    Run the GLM on several batches of data

    - `df` : the dataframe, with the data
    - `df_cache` : the dataframe, with the data
    - `cells` : dataframe describing the possible cells
    - `formula_method` : the function to create the formulae
    - `glmtool` : the tool to use
    - `Dist` : the distribution to use
    - `unitwise` : if true, we do the GLM on the unitwise data
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
        modelz=ThreadSafeDict())
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
        prog = Progress(length(glmsets); desc="GLM spike counts")
        for (indep, x_key, f) in glmsets                                                

            # Get 🔑               
            if unitwise
                unitstr = !isempty(unitrep) ?
                    replace(string(f.lhs),unitrep...) : string(f.lhs)
                unit=parse(Int, unitstr)
                out_key = (;unit, indep, x_key...)
            else
                out_key = (;indep, x_key...)
            end

            # Run model 🤖
            try
                XX, y = df_batches[x_key]
                modelz[(;indep, x_key, unit)] = glm_(f, XX, y, Dist;
                                       type=glmtool, xtrans=xtrans, ytrans=ytrans)
                cacheXY[out_key] = (;XX, y)
            catch exception
                @infiltrate
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
    function run_glm!(full_df::DataFrame, formula_method::Function,
                        cells::DataFrame, glmtool::Symbol; kws...)
        df_batches = OrderedDict(
            "batch" => (full_df, full_df)
        )   
        run_glm!(full_df, df_batches, cells, formula_method, glmtool; kws...)
    end
    """
    Runs the GLM on a grouped dataframe
    """
    function run_glm!(grouped::GroupedDataFrame, formula_method::Function,
                        cells::DataFrame, glmtool::Symbol; kws...)
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

    export uicellcols
    """
    grab columsn that correspond to cell isolated props
    """
    function uicellcols(df)
        nms = [nm for nm in names(df) if endswith(nm,"_i")]
        names_cleaned = replace.(nms, ["_i"=>""])
        nms[tryparse.(Int, names_cleaned) .!== nothing]
    end


    # STRING METHODS
    nameit(x) = lowercase(begin
        first(split(first(split(string(x),"{")), "(" ))
    end)
    linkstr(link) = link === nothing ? 
    replace(nameit(canonicallink(Dist)),"link"=>"") : link
    diststr(Dist) = Dist isa String ? Dist : nameit(Dist)
    
end
