module checkpoint

    using Serialization, DataFrames, DrWatson, Infiltrator, ArgParse
    import Arrow
    using DataStructures: OrderedDict

    import ..Timeshift
    import ..Timeshift.DIutils.Table: to_dataframe, DictOfShiftOfUnit
    import Field: ReceptiveField
    import Field.metrics: metric_ban, apply_metric_ban, unstackMetricDF
    import Load: save_table_at_path

    export ts_plotdir
    export save_mains, save_shuffles, save_fields
    export load_mains, load_shuffles, load_fields
    export mainspath, shufflespath
    export archive

    # ------------------
    # SAVING AND LOADING
    # ------------------
    """
    ts_plotdir

    verion of plotsdir that autobegins at 
    the FOLDER FOR PLOT

    plotsdir/timeshift/
    """
    function ts_plotdir(;file="", kws...)
        parent_folder = plotsdir("timeshift")
        if !(isdir(parent_folder))
            mkdir(parent_folder)
        end
        name = replace(basename(_pathshiftdat(;kws...)), ".serial"=>"")
        name = joinpath(parent_folder, name)
        if !(isdir(name))
            mkdir(name)
        end
        joinpath(name,file)
    end

    function saveTauMax_cellTable(imaxdf::DataFrame, animal::String, day::Int,
        tag::String)
        tag = "$(tag)_tauMax"
        Load.save_cell_taginfo(imaxdf, animal, day, tag)
    end

    function path(store::String; archive="", kws...)
        parent_folder = datadir("exp_pro", "timeshift")
        folder = joinpath(parent_folder, store)
        archive == "" ? folder : join([folder,archive], ".")
    end
    mainspath(;kws...)    = path("mains"; kws...)
    shufflespath(;kws...) = path("shuffles"; kws...)
    fieldspath(;kws...)   = path("fields"; kws...)
    shufflefieldspath(;kws...)   = path("shufflefields"; kws...)

    storekeys(store::String;kws...) = keys(deserialize(path(store;
                                                                 kws...)))
    mainskeys(;kws)    = storekeys("mains";kws...)
    shuffleskeys(;kws) = storekeys("shuffles";kws...)
    fieldskeys(;kws)   = storekeys("fields";kws...)

    function save_store(store, data; dataframe::Bool=false)
        if dataframe
            name = path(store) * "_dataframe"
            save_table_at_path(data, name, "arrow")
        else
            name = path(store)
            if overwrite
                serialize(name, data)
            else
                store = load_store(store)
                data  = merge(store, data)
                serialize(name, data)
            end
        end
    end

    function load_store(store; dataframe::Bool=false, archive="")
        if dataframe
            name = path(store;archive) * "_dataframe.arrow"
            println("Loading $name")
            DataFrame(Arrow.Table(name))
        else
            name = path(store;archive)
            println("Loading $name")
            D = deserialize(name)
        end
    end
    load_mains(;kws...)         = load_store("mains"; kws...)
    load_shuffles(;kws...)      = load_store("shuffles"; kws...)
    load_fields(;kws...)        = load_store("fields"; kws...)
    load_shufflefields(;kws...) = load_store("shufflefields"; kws...)

    export change_store
    """
        change_store(store, λk::Function=identity, λv::Function=identity; filt::Function=x->true, stoponfirst::Bool=true)
    
    change a store contents easily
    """
    function change_store(store, λk::Function=identity, λv::Function=identity; filt::Function=x->true, stoponfirst::Bool=true)
        store = load_store(store; dataframe=false)
        store = filter(filt, store)
        for (i,(key,value)) in enumerate(store)
            if stoponfirst && i == 1
                println("Before")
                @infiltrate
            end
            oldkey, newkey = key, λk(key)
            newval = λv(value)
            store[newkey] = newval
            if oldkey != newkey
                pop!(store,oldkey)
            end
            if stoponfirst && i == 1
                println("After")
                @infiltrate
            end
        end
        store
    end


    """
    old = OrderedCollections.OrderedDict{Float64, Dict{Any(NamedTuple), Any}} 
    new = OrderedCollections.OrderedDict{Float64, OrderedCollections.OrderedDict{NamedTuple,AdaptiveRF}}
    """
    function save_mains(M::AbstractDict; overwrite::Bool=false)
        if any(x->x === nothing, values(M))
            M = typeof(M)(k=>v for (k,v) in M if v !== nothing)
        end
        M = cut_the_fat(M)
        name = mainspath()
        if isfile(name) && !(overwrite)
            @info "Preloading existing $name"
            D = deserialize(name)
        else
            D = Dict()
        end
        D = merge(D, M)
        @info "Saving $name"
        if keytype(D) == Any
            D = Dict{NamedTuple, Any}(key=>value for (key, value) in D)
        end
        serialize(name, D)
    end
    function save_mains(M::DataFrame)
        name = mainspath() * "_dataframe.arrow"
        Arrow.write(name, M)
    end
    save_mains(M::T where T<: Timeshift.DictOfShiftOfUnit) =
                              save_mains(Timeshift.ShiftedFields(M))

    """
        cut_the_fat

    In case new strutures like (ReceptiveField and ShiftedFields and DictOfShiftOfUnit)
    exist in our dict of results, transform those into raw metrics. Mains only include
    main effect metrics, not full on receptive fields.
    """
    function cut_the_fat(M::AbstractDict)
        if valtype(M)<:ReceptiveField
            R = Dict{keytype(M), DataFrame}()
        else
            R = copy(M)
        end
        for (k,v) in M
            if v isa AbstractDict
                R[k] = cut_the_fat(M[k])
            elseif typeof(v)<:ReceptiveField
                R[k] = to_dataframe(apply_metric_ban(v.metrics); explode=false)
            elseif v isa Timeshift.ShiftedField || v isa  Timeshift.ShiftedFields
                R[k] = unstackMetricDF(to_dataframe(apply_metric_ban(v.metrics); explode=false))
            elseif v isa AbstractDataFrame && !hasproperty(v, :unit)
                R[k] = unstackMetricDF(apply_metric_ban(v))
            else
                R[k] = v
            end
        end
        R
    end

    """
    """
    function save_fields(M::AbstractDict; overwrite::Bool=false, 
                                          archive::String="")
        if any(x->x === nothing, values(M))
            M = typeof(M)(k=>v for (k,v) in M if v !== nothing)
        end
        name = fieldspath(;archive)
        if isfile(name) && !(overwrite)
            @info "Preloading existing $name"
            D = deserialize(name)
        else
            D = Dict()
        end
        D = merge(D, M)
        @info "Saving $name"
        if keytype(D) == Any
            D = Dict{NamedTuple, Any}(key=>value for (key, value) in D)
        end
        serialize(name, D)
    end
    function save_fields(F::DataFrame)
        name = fieldspath() * "_dataframe.arrow"
        Arrow.write(name, F)
    end

    function save_shufflefields(S::AbstractDict; overwrite::Bool=false,
        num_examples::Int=3)
        if any(x->x .=== nothing, values(S))
            S = typeof(S)(k=>v for (k,v) in S if v !== nothing)
        end
        name = shufflefieldspath()
        S = collect_examples(S, x->hasproperty(x, :shuffle), num_examples)
        if isfile(name) && !(overwrite)
            @info "Preloading existing $name"
            D = deserialize(name)
        else
            D = Dict()
        end
        D = merge(D, S)
        @info "Saving $name"
        serialize(name, D)
    end
    
    function collect_examples(S::AbstractDict, key_lambda::Function,
            num_examples::Int)
        @error "Not implemented!"
        valid = key_lambda.(collect(keys(S)))
        if any(valid)
            keys = [key for key in keys(S) if key_lambda(key)]
            S = typeof(S)(k=>S[k] for k in Iterators.take(keys, num_examples))
        end
        R = typeof(S){valtype(S)}{keytype(S)}()
        for (k,v) in S
            R[k] = v
        end
    end


    function save_shuffles(S::AbstractDict; overwrite::Bool=false)
        if any(x->x === nothing, values(S))
            S = typeof(S)(k=>v for (k,v) in S if v !== nothing)
        end
        name = shufflespath()
        @info "save_shuffles: cut the fat"
        S = cut_the_fat(S)
        if isfile(name) && !(overwrite)
            @info "Preloading existing $name"
            D = deserialize(name)
        else
            D = Dict()
        end
        D = merge(D, S)
        @info "Saving $name"
        serialize(name, D)
    end
    function save_shuffles(S::DataFrame)
        name = shufflespath() * "_dataframe.arrow"
        Arrow.write(name, S)
    end

    function archivepath(store::String; archive=".archive")::String
        pathfunc = eval(Symbol(store * "path"))
        path = pathfunc()
        archivepath = path * archive
    end

    function archive(store::String, keysearch::NamedTuple; archive="archive", overwrite::Bool=false)
        storepath = path(store)
        archivepath = path(store; archive)
        @info "Loading archive=$storepath and archivepath=$archivepath"
        store, archive = deserialize(storepath), 
                         overwrite ? OrderedDict() : deserialize(archivepath)
        matches = match(keys(store), keysearch)
        @infiltrate
        if !(isempty(matches))
            @info "Matches" matches
            for key in matches
                @info "Archiving" key
                push!(archive, key=> pop!(store, key))
            end
            archive = zip(keys(archive), values(archive))
            store = typeof(store)(key=>value for (key,value) in store)
            Base.rehash!(store)
            @info("Serializing archive")
            serialize(archivepath, archive)
            @info("Serializing store")
            serialize(storepath, store)
        else
            @info "No matches"
        end
        nothing
    end
    function unarchive(store::String, keysearch::NamedTuple; archive="archive")
        storepath   = path(store)
        archivepath = path(store; archive)
        store, archive = deserialize(storepath), 
                         deserialize(archivepath) 
        matches = match(keys(store), keysearch)
        @info "Matches" matches
        for key in matches
            @info "Archiving" key
            push!(store, pop!(archive, key))
        end
        @infiltrate
        @info("Serializing archive")
        serialize(archivepath, archive)
        @info("Serializing store")
        serialize(storepath, store)
        nothing

    end
    
    function Base.match(store::String, search::NamedTuple)
        loadfunc = eval(Symbol("load_" * store))
        I = loadfunc()
        match(I, search)
    end
    match(dict_w_ntkeys::AbstractDict, search::NamedTuple) =
                                        match(keys(dict_w_ntkeys), search)
    function Base.match(ntkeys::Union{Base.KeySet, Vector}, search::NamedTuple)
        ntkeys = ntkeys isa Base.KeySet ? collect(ntkeys) : ntkeys
        [ntkey for ntkey in ntkeys if match(ntkey, search)]
    end
    function Base.match(ntkey::NamedTuple, search::NamedTuple)
        onekey_keys = keys(ntkey)
        answer = true
        for key in keys(search)
            if !(key ∈ onekey_keys &&
                 getproperty(search, key) == getproperty(ntkey, key))
                answer = false
            end
        end
        answer
    end
    export argparse
    """
        argparse(args=nothing; return_parser::Bool=false)

    return command line parser for flags that control my manifold
    related analyses
    """
    function argparse(mod::Module, pos...; kws...)
        @eval mod begin
            if isdefined(mod, :USER_ARGS)
                argparse(mod.USER_ARGS, $(pos)...; $(kws)...) 
            else
                argparse($(pos)...; $(kws)...)
            end
        end
    end
    function argparse(args=nothing; return_parser::Bool=false)
        parser = ArgParseSettings()
        @add_arg_table parser begin

            # TIMESHIFTING PARAMS
            "--init_shift", "-i"
                help="Intitial explored shift"
                default=-2.0f0
                arg_type=Float32
            "--final_shift", "-f"
                help="Final explored shift"
                arg_type=Float32
                default=2.0f0
            "--period_shift", "-p"
                help="Period of shift"
                arg_type=Float32
                default=0.05f0

            # ADAPTIVE GRIDING 
            # Width of field params
            "--width", "-w"
                help="Default field spacing widths"
                default=[5.0f0]
                arg_type=Vector{Float32}
            "--thresh", "-t"
                help="Threshold of seconds of time to accumulate per adaptive bin"
                default=[1.5f0] # seconds
                arg_type=Vector{Float32}
            "--props", "-P"
                help="Which properties to use in the receptive fields"
                default=["x", "y"]
                arg_type=Vector{String}

            # SPIKE FILTErings
            "--spike_filt", "-S"
                help="filtration options for spikes"
                default=""
                arg_type=String

            # CONSTRAINTS: Some scripts operate on all of these
            "--animal", "-A"
                help = "the animal to run default: the super animal"
                arg_type = String
                default = "RY16"
            "--day", "-D"
                help = "the day to run"
                arg_type = Int
                default = 36
            "--dataset"
                # In the future, I'm intending this to point into a GoalFetch module
                 help = "dataset preset"
                 arg_type = Int
                 default = 0
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
    function postprocess(opt)
        if opt["spike_filt"] == ""
            opt["spike_filt"] = nothing
        else
            opt["spike_filt"] = Symbol(spike_filt)
        end
        opt
    end

    # ===================
    # Old Ways
    # ==================
    #"""
    #_pathshiftdat

    #provides path for saved data

    #TODO make more flexible
    #"""
    #function _pathshiftdat(;metric=nothing, shifts=nothing, 
    #            tag="", props=nothing, splitby=nothing, kws...)

    #        if props isa Nothing
    #            @error "field keyword props must be given"
    #        end
    #        if splitby isa Nothing
    #            @error "field keyword splitby must be given"
    #        end

    #        parent_folder = datadir("exp_pro", "timeshift")
    #        if !(isdir(parent_folder))
    #            mkdir(parent_folder)
    #        end

    #        tag = isempty(tag) ? tag : "_$tag"
    #        if shifts !== nothing
    #            start, stop = round(shifts[begin],digits=3),
    #                          round(shifts[end],  digits=3)
    #            N = length(shifts)
    #            shifts = "_Nstartstop=($N,$start:$stop)"
    #        else
    #            @error "No shifts name provided"
    #        end
    #        if metric === nothing
    #            @warn "No metric name provided, assuming metric='field'"
    #            metric = "field"
    #        end
    #        jf(x) = join(x,'-')
    #        props   = "props=$(jf(props))"
    #        splitby = "_splitby=$(jf(splitby))"
    #        joinpath(parent_folder, "$props$splitby$shifts$tag.serial")
    #    end
    #function saveshifts(main=nothing, shuffle=nothing; overwrite=false, kws...)

    #    name = _pathshiftdat(;kws...)
    #    if isfile(name)
    #        D = deserialize(name)
    #    else
    #        D = Dict()
    #    end

    #    kws = Dict(zip(keys(kws), values(kws)))

    #    # Remove parts that will break serialization
    #    if :filters ∈ propertynames(kws)
    #        pop!(kws, :filters)
    #    end
    #    if :metric ∈ kws
    #        metric = pop!(kws, :metric)
    #    else
    #        metric = nothing
    #    end
    #    if :shifts ∈ kws
    #        shifts = pop!(kws, :shifts)
    #    else
    #        @error "Provide shifts"
    #    end

    #    d = Dict(:main     => main,
    #             :shuffle  => shuffle,
    #             :shifts   => shifts,
    #             :metric   => metric,
    #             :fieldkws => kws)
    #    D = overwrite ? merge(D, d) : merge(d, D)
    #    serialize(name, D)
    #end
    #function loadshifts(;kws...)::Dict
    #    name = _pathshiftdat(;kws...)
    #    println(name)
    #    deserialize(name)
    #end

end
