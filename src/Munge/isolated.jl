module isolated
    using JLD2, ArgParse, DrWatson, DIutils.dict, RecipesBase, GLM, 
          DataFramesMeta, Statistics, NaNStatistics, Infiltrator
    using DataStructures: OrderedDict
    import DIutils: Table

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
        load_iso

    Returns a dictionary of the isolated variables of interest
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

    export construct_predict_spikecount
    function construct_predict_spikecount(df, cells, independent_area="CA1";
            other_vars=[], other_ind_vars=[])
        uArea = unique(cells.area)
        @assert length(uArea) == 2 "Only supports two area dataframes"
        dependent_area = setdiff(uArea, [independent_area])

        dep_neurons = @subset(cells,:area .==dependent_area).unit
        ind_neurons = @subset(cells,:area .==independent_area).unit
        filter!(n->string(n) ∈ names(df), dep_neurons)
        filter!(n->string(n) ∈ names(df), ind_neurons)
        
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

    export construct_predict_iso
    function construct_predict_iso(df, cells, independent_area="CA1", type=:has;
            other_vars=[], other_ind_vars=[])
        uArea = unique(cells.area)
        @assert length(uArea) == 2 "Only supports two area dataframes"
        ind_neurons = @subset(cells,:area .==independent_area).unit
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

    """
    Table.to_dataframe
    """
    function Table.to_dataframe(D::AbstractDict, func::Function)::DataFrame
        kt = keytype(D)
        D = Dict{kt, Any}(k=>func(k, v) for (k,v) in D)
        Table.to_dataframe(D)
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
        out
    end
end
