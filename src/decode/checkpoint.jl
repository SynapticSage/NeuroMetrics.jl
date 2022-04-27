
using HDF5
using Arrow

function save_checkpoint(Module m)
    cp(x, type="arrow") = joinpath(dirname(decode_file), "$(x).$type")
    @info "Checkpointing decode variables"
    begin
        Arrow.write(cp("cells"), m.cells)
        Arrow.write(cp("spikes"), m.spikes)
        Arrow.write(cp("beh"), m.beh)
        Arrow.write(cp("lfp"), m.lfp)
        Arrow.write(cp("cycles"), m.cycles)
        Arrow.write(cp("ripples"), m.ripples)
        h5open(cp("decode","h5"), "w") do file
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

function load_checkpoint()
    @error "Have not written for load yet"
    begin
        Arrow.write(cp("cells"), m.cells)
        Arrow.write(cp("spikes"), m.spikes)
        Arrow.write(cp("beh"), m.beh)
        Arrow.write(cp("lfp"), m.lfp)
        Arrow.write(cp("cycles"), m.cycles)
        Arrow.write(cp("ripples"), m.ripples)
        h5open(cp("decode","h5"), "w") do file
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
