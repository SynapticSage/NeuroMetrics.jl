module checkpoint

    using Serialization
    using DataFrames
    using DrWatson
    import Arrow

    export ts_plotdir
    export save_mains, save_shuffles, save_fields
    export load_mains, load_shuffles, load_fields
    export mainspath, shufflespath

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

    function saveTauMax_cellTable(df_imax::DataFrame, animal::String, day::Int,
        tag::String)
        tag = "$(tag)_tauMax"
        Load.save_cell_taginfo(df_imax, animal, day, tag)
    end

    function mainspath()
        parent_folder = datadir("exp_pro", "timeshift")
        name = joinpath(parent_folder, "mains")
    end

    function save_mains(M::AbstractDict)
        name = mainspath()
        if isfile(name)
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

    function save_fields(M::AbstractDict)
        name = fieldspath()
        if isfile(name)
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

    function load_mains(;dataframe::Bool=false)
        if dataframe
            name = mainspath() * "_dataframe.arrow"
            DataFrame(Arrow.Table(name))
        else
            name = mainspath()
            D = deserialize(name)
        end
    end

    function fieldspath()
        parent_folder = datadir("exp_pro", "timeshift")
        name = joinpath(parent_folder, "fields")
    end

    function load_fields(;dataframe::Bool=false)
        if dataframe
            name = fieldspath() * "_dataframe.arrow"
            DataFrame(Arrow.Table(name))
        else
            name = fieldspath()
            D = deserialize(name)
        end
    end


    function shufflespath()
        parent_folder = datadir("exp_pro", "timeshift")
        name = joinpath(parent_folder, "shuffles")
    end

    function save_shuffles(S::AbstractDict)
        name = shufflespath()
        if isfile(name)
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

    function load_shuffles(;dataframe::Bool=false)
        if dataframe
            name = shufflespath() * "_dataframe.arrow"
            DataFrame(Arrow.Table(name))
        else
            name = shufflespath()
            D = deserialize(name)
        end
    end

    # ===================
    # Old Ways
    # ==================
    """
    _pathshiftdat

    provides path for saved data

    TODO make more flexible
    """
    function _pathshiftdat(;metric=nothing, shifts=nothing, 
                tag="", props=nothing, splitby=nothing, kws...)

            if props isa Nothing
                @error "field keyword props must be given"
            end
            if splitby isa Nothing
                @error "field keyword splitby must be given"
            end

            parent_folder = datadir("exp_pro", "timeshift")
            if !(isdir(parent_folder))
                mkdir(parent_folder)
            end

            tag = isempty(tag) ? tag : "_$tag"
            if shifts != nothing
                start, stop = round(shifts[begin],digits=3),
                              round(shifts[end],  digits=3)
                N = length(shifts)
                shifts = "_Nstartstop=($N,$start:$stop)"
            else
                @error "No shifts name provided"
            end
            if metric == nothing
                @warn "No metric name provided, assuming metric='field'"
                metric = "field"
            end
            jf(x) = join(x,'-')
            props   = "props=$(jf(props))"
            splitby = "_splitby=$(jf(splitby))"
            name = joinpath(parent_folder, "$props$splitby$shifts$tag.serial")
        end
    function saveshifts(main=nothing, shuffle=nothing; overwrite=false, kws...)

        name = _pathshiftdat(;kws...)
        if isfile(name)
            D = deserialize(name)
        else
            D = Dict()
        end

        kws = Dict(zip(keys(kws), values(kws)))

        # Remove parts that will break serialization
        if :filters ∈ propertynames(kws)
            pop!(kws, :filters)
        end
        if :metric ∈ kws
            metric = pop!(kws, :metric)
        else
            metric = nothing
        end
        if :shifts ∈ kws
            shifts = pop!(kws, :shifts)
        else
            @error "Provide shifts"
        end

        d = Dict(:main     => main,
                 :shuffle  => shuffle,
                 :shifts   => shifts,
                 :metric   => metric,
                 :fieldkws => kws)
        D = overwrite ? merge(D, d) : merge(d, D)
        serialize(name, D)
    end
    function loadshifts(;kws...)::Dict
        name = _pathshiftdat(;kws...)
        println(name)
        D = deserialize(name)
    end


end
