module checkpoint

    using HDF5
    using Arrow
    export HDF5
    export Arrow
    using ProgressMeter
    using Infiltrator
    using Glob
    using DataFrames

    pathname(x, decode_file, type="arrow") = joinpath(dirname(decode_file), "split=$(extract_splitnum(decode_file))_$(x).$type")

    export load_checkpoints, save_checkpoint, load_checkpoint, extract_splitnum

    function save_checkpoint(m::Module, decode_file; split, overwrite=true)
        @info "Checkpointing decode variables"
        pn(x)  = pathname(x, decode_file)
        h5file = pathname("split=$(split)_decode",decode_file, "h5")
        if overwrite == false && isfile(h5file)
            @info "File exists and returning"
            return nothing
        end
        begin
            Arrow.write(pn("split=$(split)_cells"),   m.cells)
            Arrow.write(pn("split=$(split)_spikes"),  m.spikes)
            Arrow.write(pn("split=$(split)_beh"),     m.beh)
            Arrow.write(pn("split=$(split)_lfp"),     m.lfp)
            Arrow.write(pn("split=$(split)_cycles"),  m.cycles)
            Arrow.write(pn("split=$(split)_ripples"), m.ripples)
            @info "Saving $h5file"
            h5open(h5file, "w") do file
                create_group(file, "decode")
                println(file["decode"])
                file["decode/dat"] = m.dat
                #if hasproperty(m, :ripple)
                #    file["decode/ripple"] = m.ripple
                #    file["decode/theta"] = m.theta
                #    file["decode/non"] = m.non
                #end
                file["decode/x"]   = m.x
                file["decode/y"]   = m.y
                file["decode/T"]   = m.T
            end
            nothing
        end
    end


    function get_arrow_vars()
        [:cells, :spikes, :beh, :lfp, :cycles, :ripples]
    end
    function get_h5_vars()
        [:dat, :theta, :ripple, :non, :x, :y, :T]
    end

    function extract_splitnum(file::String)
        file = basename(file)
        file = split(file, ".")[begin]
        parse(Int, split(split(file, "split=")[end],"_")[begin])
    end

    function load_checkpoint(decode_file::String; vars=nothing)
        pn(x) = pathname(x, decode_file) 
        D = Dict()
        arrow_vars  = get_arrow_vars()
        h5_vars     = get_h5_vars()
        if vars !== nothing
            arrow_vars = filter(x->x∈vars, arrow_vars)
            h5_vars    = filter(x->x∈vars, h5_vars)
        end

        for var in arrow_vars
            D[var]    = DataFrame(Arrow.Table(pn(String(var))))
        end
        filename = if isfile(decode_file)
            decode_file
        elseif isfile(pathname("decode", decode_file, "h5"))
            pathname("decode", decode_file, "h5")
        end
        h5open(filename, "r") do file
            for var in h5_vars
                if "$(String(var))" in keys(file["decode"])
                    D[var]    = read(file, "decode/$(String(var))")
                end
            end
        end
        return D
    end

    """
    Load into another module

    HOw woulld I do this without globals?
    """
    function load_checkpoint(m::Module, decode_file::String; vars=nothing)
        D = load_checkpoint(decode_file; vars=vars)
        global D
        @eval m global D
        @eval m global key
        for (key,value) in D
            @eval m Meta.parse("$key = D[:$key]")
        end
        nothing
    end

    """
    variables not to be concatonated when multiple checkpoint files
    """
    noncatvars = [:x, :y]


    function generalize_path(decode_file)
        dir = dirname(decode_file)
        base = split(basename(decode_file),".")[1]
        base = split(base, "=")[1]
        path = joinpath(dir, base)
    end

    """
    load_checkpoints

    way to load multiple checkpoints concatonated along an axis
    """
    function load_checkpoints(decode_file::String;
            vars::Union{Nothing, Vector{Symbol}}=nothing)

        path = generalize_path(decode_file)
        dir  = dirname(path)
        base = basename(path)
        @debug "dir=$dir, base=$base"

        D = Dict()
        L = Vector{Dict}()
        @showprogress for (i,file) in enumerate(glob(base.*"=*.nc", dir))
            @debug "file=$file"
            try
                @assert length(split(basename(file), "split=")) < 3
            catch
                @infiltrate
            end
            tmp = load_checkpoint(file, vars=vars)
            tmp = Load.fix_complex(tmp)
            tmp = Load.fix_rgba(tmp)
            if i == 1
                D = merge(D,tmp)
            else
                push!(L, tmp)
            end
            tmp = nothing
            
            GC.gc(true)
        end
        for field in keys(D)
            if field ∈ noncatvars
            elseif field ∈ get_arrow_vars() || ndims(D[field]) == 1
                D[field] = cat(D[field], [l[field] for l in L]..., dims=1)
            elseif field ∈ get_h5_vars()
                D[field] = cat(D[field], [l[field] for l in L]..., dims=3)
            end
            [delete!(l,field) for l in L]
            nothing
        end

        D
    end

end
