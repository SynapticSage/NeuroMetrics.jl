
using HDF5
using Arrow
export HDF5
export Arrow
using Glob

pathname(x, decode_file, type="arrow") = joinpath(dirname(decode_file), "split=$(extract_splitnum(decode_file))_$(x).$type")
export load_checkpoints

function save_checkpoint(m::Module, decode_file; split, overwrite=true)
    @info "Checkpointing decode variables"
    pn(x) = pathname(x, decode_file)
    h5file = pathname("split=$(split)_decode",decode_file, "h5")
    if overwrite == false && isfile(h5file)
        @info "File exists and returning"
        return nothing
    end
    begin
        Arrow.write(pn( "split=$(split)_cells"),   m.cells)
        Arrow.write(pn( "split=$(split)_spikes"),  m.spikes)
        Arrow.write(pn( "split=$(split)_beh"),     m.beh)
        Arrow.write(pn( "split=$(split)_lfp"),     m.lfp)
        Arrow.write(pn( "split=$(split)_cycles"),  m.cycles)
        Arrow.write(pn( "split=$(split)_ripples"), m.ripples)
        @info "Saving $h5file"
        h5open(h5file, "w") do file
            create_group(file, "decode")
            println(file["decode"])
            file["decode/dat"] = m.dat
            file["decode/ripple"] = m.ripple
            file["decode/theta"] = m.theta
            file["decode/non"] = m.non
            file["decode/x"]   = m.x
            file["decode/y"]   = m.y
            file["decode/T"]   = m.T
        end
        nothing
    end
end


function get_arrow_vars()
    arrow_vars  = [:cells, :spikes, :beh, :lfp, :cycles, :ripples]
end
function get_h5_vars()
    h5_vars     = [:dat, :theta, :ripple, :non, :x, :y, :T]
end

function extract_splitnum(file::String)
    file = basename(file)
    file = split(file, ".")[begin]
    num = parse(Int, split(file, "split=")[end])
end


function load_checkpoint(decode_file::String; vars=nothing)
    @infiltrate
    pn(x) = pathname(x, decode_file) 
    D = Dict()
    arrow_vars  = get_arrow_vars()
    h5_vars     = get_h5_vars()
    if vars != nothing
        arrow_vars = filter(x->x∈vars, arrow_vars)
        h5_vars    = filter(x->x∈vars, h5_vars)
    end

    for var in arrow_vars
        D[var]    = DataFrame(Arrow.Table(pn(String(var))))
    end
    h5open(pathname("decode", decode_file, "h5"), "r") do file
        for var in h5_vars
            D[var]    = read(file, "decode/$(String(var))")
        end
    end
    return D
end

function load_checkpoint(m::Module, decode_file::String; vars=nothing)
    D = load_checkpoint(decode_file; vars=vars)
    @error "Not implemented here"
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
function load_checkpoints(decode_file::String; vars::Union{Nothing, Vector{Symbol}}=nothing)

    path = generalize_path(decode_file)
    dir  = dirname(path)
    base = basename(path)
    @debug "dir=$dir, base=$base"

    for (i,file) in enumerate(glob(base.*"=*.nc", dir))
        @debug "file=$file"
        tmp = load_checkpoint(file, vars=vars)
        if i == 1
            D = tmp
        else
            for field in keys(tmp)
                if field ∈ noncatvars
                    D[x] = tmp[x]
                elseif field ∈ get_arrow_vars() || ndims(tmp[x]) == 1
                    D[x] = cat(D[x], tmp[x], dims=1)
                elseif field ∈ get_h5_vars()
                    D[x] = cat(D[x], tmp[x], dims=3)
                end
            end
        end
        D
    end

end
