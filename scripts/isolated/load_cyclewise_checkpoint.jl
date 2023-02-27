using DrWatson; quickactivate("goal-code")
include(scriptsdir("isolated","imports_isolated.jl"))
include("./imports_isolated.jl")
opt = parser()
opt["cycles"] = 8
# Data
@time data = load_iso(opt)
lfp       = data["lfp"]; @assert(lfp isa DataFrame)
spikes    = data["spikes"] ; @assert(spikes isa DataFrame)
allspikes = data["allspikes"]; @assert(allspikes isa DataFrame)
tsk       = data["tsk"]; @assert(tsk isa DataFrame)
cells     = data["cells"]; @assert(cells isa DataFrame)
beh       = data["beh"]; @assert(beh isa DataFrame)
cycles    = data["cycles"]; @assert(cycles isa DataFrame)
Rdf       = data["Rdf"]; @assert(Rdf isa DataFrame)
beh.speed = abs.(beh.smoothvel)

fn = path_iso(opt; append="_cyclewise")
grd = occ = df = nothing
jldopen(fn, "r") do storage
    @eval Main grd = $(storage["grd"])
    @eval Main occ = $(storage["occ"])
    @eval Main df  = $(storage["df"])
end

function get_dx_dy(df, relcyc)
    dx = @subset(df, :relcycs .== relcyc)
    dy = @subset(df, :relcycs .== 0)
    dx = groupby(dx, [:cyc_batch, :cyc_match])
    dy = groupby(dy, [:cyc_batch, :cyc_match])
    kx, ky = keys(dx.keymap), keys(dy.keymap)
    k = intersect(kx,ky)
    dx = sort(vcat([dx[kk] for kk in k]...), [:cyc_batch, :cyc_match])
    dy = sort(vcat([dy[kk] for kk in k]...), [:cyc_batch, :cyc_match])
    dx, dy
end
@time dx_dy = Dict(relcyc => get_dx_dy(df, relcyc) for relcyc 
    in -opt["cycles"]:opt["cycles"])
DIutils.pushover("Finished load_cyclewise_checkpoint")
