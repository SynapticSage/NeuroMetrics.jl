include("imports_isolated.jl")
# begin

    # Parse the command line
    opt = parser()
    # Data
    data = load_iso(opt)
    lfp       = data["lfp"];     @assert(lfp isa DataFrame)
    spikes    = data["spikes"] ; @assert(spikes isa DataFrame)
    allspikes = data["allspikes"]; @assert(allspikes isa DataFrame)
    tsk       = data["tsk"];     @assert(tsk isa DataFrame)
    cells     = data["cells"];   @assert(cells isa DataFrame)
    beh       = data["beh"];     @assert(beh isa DataFrame)
    cycles    = data["cycles"];  @assert(cycles isa DataFrame)
    Rdf       = data["Rdf"];     @assert(Rdf isa DataFrame)
    beh.speed = abs.(beh.smoothvel)

# end
save_kws = (;pfc_rate_analy=true)
filt = Filt.get_filters()

datacut = opt["filt"]
animal  = opt["animal"]
day     = opt["day"]
tet     = opt["tet"]
Plot.setappend((;animal, day, tet))
Munge.nonlocal.setunfilteredbeh(DI.load_behavior(animal,day);
                               animal, day)


# CONSTANTS and trackers
# --------
isonames  = OrderedDict(false => :adjacent, true=>:isolated)
filt_desc = OrderedDict(:all => "> 2cm/s")
clab      = OrderedDict(-1 => "nontask", 
                         0 => "cue", 
                         1=> "mem", 
                         missing=>"sleep")
val = :value
opt["matched"] = 3    # how many matched non-iso cycles to add per iso cycle
opt["has_df"] = false # tracks state about whether df var for glm is backed up
                      # in a jld2 file
matchprops = [:x, :y, :speed, :startWell, :stopWell]
Munge.nonlocal.setclab(clab)
