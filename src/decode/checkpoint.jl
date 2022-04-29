
using HDF5
using Arrow
export HDF5
export Arrow

pathname(x, decode_file, type="arrow") = joinpath(dirname(decode_file), "$(x).$type")

function save_checkpoint(m::Module, decode_file)
    @info "Checkpointing decode variables"
    pn(x) = pathname(x, decode_file)
    begin
        Arrow.write(pn("cells"),   m.cells)
        Arrow.write(pn("spikes"),  m.spikes)
        Arrow.write(pn("beh"),     m.beh)
        Arrow.write(pn("lfp"),     m.lfp)
        Arrow.write(pn("cycles"),  m.cycles)
        Arrow.write(pn("ripples"), m.ripples)
        h5open(pathname("decode","h5"), "w") do file
            create_group(file, "decode")
            println(file["decode"])
            file["decode/dat"] = m.dat
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
        D[:dat] = read(file, "decode/dat")
        D[:x]   = read(file, "decode/x")
        D[:y]   = read(file, "decode/y") 
        D[:T]   = read(file, "decode/T")
    end
    return D
end
