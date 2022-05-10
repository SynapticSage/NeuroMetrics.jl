
using HDF5
using Arrow
export HDF5
export Arrow

pathname(x, decode_file, type="arrow") = joinpath(dirname(decode_file), "$(x).$type")

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


function load_checkpoint(m::Module, decode_file)
    pn(x) = pathname(x, decode_file) 
    D = Dict()
    D[:cells]    = DataFrame(Arrow.Table(pn("cells")))
    D[:spikes]   = DataFrame(Arrow.Table(pn("spikes")))
    D[:beh]      = DataFrame(Arrow.Table(pn("beh")))
    D[:lfp]      = DataFrame(Arrow.Table(pn("lfp")))
    D[:cycles]   = DataFrame(Arrow.Table(pn("cycles")))
    D[:ripples]  = DataFrame(Arrow.Table(pn("ripples")))
    h5open(pathname("decode", decode_file, "h5"), "r") do file
        D[:dat]    = read(file, "decode/dat")
        D[:theta]  = read(file, "decode/theta")
        D[:ripple] = read(file, "decode/ripple")
        D[:non] = read(file, "decode/non")
        D[:x]   = read(file, "decode/x")
        D[:y]   = read(file, "decode/y") 
        D[:T]   = read(file, "decode/T")
    end
    return D
end

"""
variables not to be concatonated when multiple checkpoint files
"""
noncatvars = [:x, :y]

"""
load_checkpoints

way to load multiple checkpoints concatonated along an axis
"""
function load_checkpoints(m::Module, decode_file, vars=nothing)
end
