
    # ------------------
    # SAVING AND LOADING
    # ------------------

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

    function saveshifts(main=nothing, shuffle=nothing; overwrite=false, kws...)

        name = _pathshiftdat(;kws...)
        if isfile(name)
            D = deserialize(name)
        else
            D = Dict()
        end

        d = Dict(:main     => main,
                 :shuffle  => shuffle,
                 :shifts   => kws.shifts,
                 :metric   => kws.metric,
                 :fieldkws => fieldkws)
        D = overwrite ? merge(D, d) : merge(d, D)
        serialize(name, D)
    end
    function loadshifts(;kws...)::Dict
        name = _pathshiftdat(;kws...)
        println(name)
        D = deserialize(name)
    end

    function mainspath()
        parent_folder = datadir("exp_pro", "timeshift")
        name = joinpath(parent_folder, "mains")
    end
    function save_mains(M::AbstractDict)
        name = mainspath()
        if isfile(name)
            D = deserialize(name)
        else
            D = Dict()
        end
        D = merge(D, M)

        serialize(name, D)
    end


    function shufflespath()
        parent_folder = datadir("exp_pro", "timeshift")
        name = joinpath(parent_folder, "shuffles")
    end
    function save_shuffles(S::AbstractDict)
        name = shufflespath()
        @info "saving $name"
        if isfile(name)
            D = deserialize(name)
        else
            D = Dict()
        end
        D = merge(D, S)
        serialize(name, D)
    end
    function load_shuffles()
        name = shufflespath()
        D = deserialize(name)
    end



