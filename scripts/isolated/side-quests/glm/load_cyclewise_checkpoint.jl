using DrWatson; quickactivate("goal-code")
include(scriptsdir("isolated","imports_isolated.jl"))
include("./imports_isolated.jl")
opt = parser()
opt["cycles"]  = 8
opt["matched"] = 3    # how many matched non-iso cycles to add per iso cycle
opt["commits"] = false # whether commit() statements actually save vars

# Data
@time data = load_iso(opt)
lfp       = data["lfp"];       @assert(lfp isa DataFrame)
spikes    = data["spikes"] ;   @assert(spikes isa DataFrame)
allspikes = data["allspikes"]; @assert(allspikes isa DataFrame)
tsk       = data["tsk"];       @assert(tsk isa DataFrame)
cells     = data["cells"];     @assert(cells isa DataFrame)
beh       = data["beh"];       @assert(beh isa DataFrame)
cycles    = data["cycles"];    @assert(cycles isa DataFrame)
Rdf       = data["Rdf"];       @assert(Rdf isa DataFrame)
R         = data["R"];       
beh.speed = abs.(beh.smoothvel)

# 
fn = path_iso(opt; append="_cyclewise")
grd = occ = df = nothing
jldopen(fn, "r") do storage
    @eval Main grd = $(storage["grd"])
    @eval Main occ = $(storage["occ"])
    @eval Main df  = $(storage["df"])
end

dx_dy = Dict(relcyc => get_dx_dy(df, relcyc)
    for relcyc in minimum(df.relcycs):maximum(df.relcycs))

DIutils.pushover("Finished load_cyclewise_checkpoint")

# Models
model_spikecount   = initorget("model_spikecount")
model_cellcountiso = initorget("model_cellcountiso")
model_cellhasiso   = initorget("model_cellhasiso")

# Shuffles
shuffle_cellhasiso   = initorget("shuffle_cellhasiso")

# DF
df     = initorget("df"; obj=DataFrame())
df_iso = initorget("df_iso"; obj=DataFrame())

matchprops = [:x, :y, :speed, :startWell, :stopWell]
val = :value
