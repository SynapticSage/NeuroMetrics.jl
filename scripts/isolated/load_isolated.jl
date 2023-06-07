if isdefined(Main,:DrWatson)
    @eval Main include(scriptsdir("isolated","imports_isolated.jl"))
    @eval Main include(scriptsdir("isolated","setup_checkpoint.jl"))
else
    @eval Main include("imports_isolated.jl")
    @eval Main include("setup_checkpoint.jl")
end

# # begin
#
#     # Parse the command line
#     opt = parser()
#     # Data
#     data = load_iso(opt)
#     lfp       = data["lfp"];       @assert(lfp isa DataFrame)
#     spikes    = data["spikes"] ;   @assert(spikes isa DataFrame)
#     allspikes = data["allspikes"]; @assert(allspikes isa DataFrame)
#     tsk       = data["tsk"];       @assert(tsk isa DataFrame)
#     cells     = data["cells"];     @assert(cells isa DataFrame)
#     beh       = data["beh"];       @assert(beh isa DataFrame)
#     cycles    = data["cycles"];    @assert(cycles isa DataFrame)
#     Rdf       = data["Rdf"];       @assert(Rdf isa DataFrame)
#     beh.speed = abs.(beh.smoothvel)
#
# # end
# save_kws = (;pfc_rate_analy=true)
# filt = Filt.get_filters()

datacut = opt["filt"]
animal  = opt["animal"]
day     = opt["day"]
tetrode_set = opt["tetrode_set"]
opt["matched"] = 3    # how many matched non-iso cycles to add per iso cycle
opt["has_df"] = false # tracks state about whether df var for glm is backed up
                      # in a jld2 file
opt["commits"] = false # whether commit() statements actually save vars
val = :value # TODO these should be in opt
matchprops = [:x, :y, :speed, :startWell, :stopWell] # TODO
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
neuron_names = OrderedDict(false=>"pyr", true=>"int")
Munge.nonlocal.setclab(clab)

println("Imports isolated: ", opt)
Plot.setappend((;animal, day, tetrode_set=opt["tetrode_set"]))
DIutils.pushover("Imports isolated: $opt")

# spikes = subset(SPIKES, :animal=>a->a.==animal, :day=>d->d.==day)
# cells  = subset(CELLS, :animal=>a->a.==animal, :day=>d->d.==day)
# spikes = SPIKES
# cells = CELLS
