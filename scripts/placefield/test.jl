#using DaWatson
quickactivate(expanduser("~/Projects/goal-code"))
# Grab our raw data
using DataFrames
using KernelDensity, Distributions
using Plots, Measures
using ProgressMeter
includet(srcdir("raw.jl"))
includet(srcdir("field.jl"))
includet(srcdir("filt.jl"))
spikes, beh, ripples, cells = raw.load("RY16", 36);
@assert ("x" ∈ names(beh)) "Fuck"
function sf(p, loc)
    p.attr[:size] = (1900, 1900)
    savefig(p, loc)
    savefig(p*".svg", loc)
end
resolution = 60; # field resolution
splitby=["unit", "area"]
kws=(;resolution, splitby, filters=merge(filt.speed, filt.cellcount))

#                                                   
# ,---.|                       ,---.o     |        |
# |---'|    ,---.,---.,---.    |__. .,---.|    ,---|
# |    |    ,---||    |---'    |    ||---'|    |   |
# `    `---'`---^`---'`---'    `    ``---'`---'`---'
newkws = (; kws..., resolution=2*40, gaussian=2.3*0.5,
          filters=merge(filt.speed_lib, filt.cellcount))
place = field.get_fields(beh, spikes; newkws...);

# Generate all of the fields >  5cm/s
# --------------------------
a=field.plot.show_fields(place.hist)
b=field.plot.show_fields(place.kde)
grouping = field.group([place.hist,place.kde],["hist","kd"])
overall = field.plot.show_fieldgroups(grouping)
grouping = field.group([place.hist,place.kde],["hist","kd"], as=Dict)
name=plotsdir("fields", "individual_pf_vel=gt4")
mkdir(name)
name=plotsdir("fields", "individual_pf_vel=gt4", "place")
individualgroups = field.plot.show_fieldgroups(grouping, as=Dict,
                                              nplots_factor=0.4)
field.plot.save_dict_fieldplots(individualgroups, name)

## ,---.          |        ,---.o     |        |
# |  _.,---.,---.|        |__. .,---.|    ,---|
# |   ||   |,---||        |    ||---'|    |   |
# `---'`---'`---^`---'    `    ``---'`---'`---'
goalprops = ["currentPathLength","currentAngle"];
filters = merge(filt.speed_lib, 
                filt.notnan("currentPathLength"), 
                filt.cellcount,
                filt.correct,
                filt.notnan("currentAngle"), 
                filt.minmax("currentPathLength", 4, 90))
newkws = (;kws..., resolution=2*40, gaussian=3*0.5,
          filters=filters)
goalpath = field.get_fields(beh, spikes; props=goalprops, newkws...)
a=field.plot.show_fields(goalpath.hist)
b=field.plot.show_fields(goalpath.kde)
grouping = field.group([goalpath.hist,goalpath.kde],["hist","kd"])
overall = field.plot.show_fieldgroups(grouping)
name=plotsdir("fields", "individual_gf_vel=gt2")
mkdir(name)
name=plotsdir("fields", "individual_gf_vel=gt2", "goal")
individualgroups = field.plot.show_fieldgroups(grouping, as=Dict,
                                               nplots_factor=0.4);
field.plot.save_dict_fieldplots(individualgroups, name)
sp = raw.filterTables(spikes, filters=filters)[1]
summary = combine(groupby(sp,:unit), nrow)
sum(summary.nrow .< 100)

#goalprops = ["currentEucDist","currentAngle"];
#newkws = (;kws..., resolution=40, 
#          filters=merge(filt.speed, filt.notnan("currentEucDist"),
#                        filt.notnan("currentAngle"),
#                        filt.minmax("currentEucDist", 0, 100) ))
#fields, kdfields, xy = field.get_fields(beh, spikes; props=goalprops, newkws...)
#a=field.plot.show_fields(fields);
#b=field.plot.show_fields(kdfields);
#grouping = field.group([fields,kdfields],["hist","kd"], requireAll=true)
#overall = field.plot.show_fieldgroups(grouping)
#goaleuckd, goaleuchist = kdfields, fields;

#                                                            
# ,---.     |    o|        |                            |    
# `---.,---.|    .|---     |---.,   .    ,---.,---.,---.|    
#     ||   ||    ||        |   ||   |    |   ||   |,---||    
# `---'|---'`---'``---'    `---'`---|    `---|`---'`---^`---'
#      |                        `---'    `---'               
goalprops = ["currentPathLength","currentAngle"];
newkws = (;kws..., splitby=["stopWell", "unit","area"],
          resolution=2*40, gaussian=3*0.5,
          filters=merge(filt.speed_lib, Dict("stopWell"=>x->x.>0),
                        filt.notnan("currentPathLength"),
                        filt.correct,
                        filt.cellcount,
                        filt.notnan("currentAngle"),
                        filt.minmax("currentPathLength", 4, 90)))
goalpath = field.get_fields(beh, spikes; props=goalprops, newkws...)
#a = field.plot.show_fields(goalpath.hist)
#b = field.plot.show_fields(goalpath.kde)
grouping = field.group([goalpath.hist,goalpath.kde],["hist","kd"])
#overall = field.plot.show_fieldgroups(grouping)
name=plotsdir("fields", "individual_gfSplitBy_vel=gt2")
mkdir(name)
name=plotsdir("fields", "individual_gfSplitBy_vel=gt2", "goal")
individualgroups = field.plot.show_fieldgroups(grouping, as=Dict,
                                               nplots_factor=0.4);
field.plot.save_dict_fieldplots(individualgroups, name, ext=["svg","png","pdf"])


#                                                       
# ,---.          |    |             |                   
# |  _.,---.,---.|    |        ,---.|    ,---.,---.,---.
# |   ||   |,---||    |        |   ||    ,---||    |---'
# `---'`---'`---^`---'`---'    |---'`---'`---^`---'`---'
#                              |                        
grouping = field.group([place.hist,goalpath.hist, place.hist, goalpath.kde],
                       ["place","goal","placekde","goalkde"],
                       requireAll=true)
overall = field.plot.show_fieldgroups(grouping)
name=plotsdir("fields", "individual_gfpf_vel=gt2")
mkdir(name)
name=plotsdir("fields", "individual_gfpf_vel=gt2", "goal")
individualgroups = field.plot.show_fieldgroups(grouping, as=Dict,
                                               nplots_factor=0.4);
field.plot.save_dict_fieldplots(individualgroups, name, ext=["svg","png","pdf"])
#grouping = field.group([placehist,goalhist],["place","goal"], requireAll=true)
#overallhist = field.plot.show_fieldgroups(grouping)

# GOAL / PLACE
#goalprops = ["x","y","currentAngle"];
#newkws = (;kws..., resolution=[40, 40, 6], 
#          filters=merge(filt.speed_lib, filt.notnan("currentPathLength"),
#                        filt.notnan("currentAngle"), filt.minmax("currentPathLength", 0, 250)))
#goalplace = field.get_fields(beh, spikes; props=goalprops, newkws...)
#place = field.get_fields(beh, spikes; props=goalprops[1:2], newkws...)
#fields = field.operation.binary(goalplace.hist, place.hist)
#
#variation(x) = median(abs.(x .- median(x)))
#bootvariation(x) = statistic.boot(collect(skipnan(x)), 10000, variation)
#function bootvariation_3d(field)
#    fieldnew = zeros(size(field,1), size(field,2), 3)
#    @showprogress for i in 1:size(field,1)
#        for j in 1:size(field,2)
#            if !(all(isnan.(field[i,j,:])))
#                fieldnew[i,j,:] .= bootvariation(field[i,j,:])[1]
#            end
#        end
#    end
#end
#field_dists = field.operation.unary(fields, bootvariation_3d)

# ---------------------------------------------------------------------------
# CUE VERS MEMORY
# VEL > 1

memfields, memkdfields, grid = field.get_fields(beh, spikes; filters=filters)
cuefields, cuekdfields, grid = field.get_fields(beh, spikes; filters=filters)

# Generate all of the fields
# --------------------------
# Let's detect what's there
a=field.show_all_fields(cuekdfields, cells; textcolor=:black)
b=field.show_all_fields(memkdfields, cells; textcolor=:black)
overall = Plots.plot(a,b)
plots = plotsdir.("fields", "goalfield" .* ["cue", "mem", "cuemem"] .* "kde_vel>1")
sf(a, plots[1])
sf(b, plots[2])
sf(overall, plots[3])
a=field.show_all_fields(cuefields, cells; textcolor=:black)
b=field.show_all_fields(memfields, cells; textcolor=:black)
overall = Plots.plot(a,b)
plots = plotsdir.("fields", "goalfield" .* ["cue", "mem", "cuemem"] .* "hist_vel>1")
sf(a, plots[1])
sf(b, plots[2])
sf(overall, plots[3])

#                                              
# ,---.          |        ,---.o     |        |
# |  _.,---.,---.|        |__. .,---.|    ,---|
# |   ||   |,---||        |    ||---'|    |   |
# `---'`---'`---^`---'    `    ``---'`---'`---'
# ANY SPEED

filters = Dict("velVec"=> x->abs.(x).>=0,
               "cuemem"=> x->x.==1)
memfields, memkdfields = field.get_fields(beh, spikes;
                                               filters=filters,
                                               props=goalprops
                                              )
filters = Dict("velVec"=> x->abs.(x).>=0,
               "cuemem"=> x->x.==0)
cuefields, cuekdfields = field.get_fields(beh, spikes; 
                                               filters=filters,
                                               props=goalprops
                                              )

# Generate all of the fields
# --------------------------
# Let's detect what's there
a=field.show_all_fields(cuekdfields, cells; textcolor=:black)
b=field.show_all_fields(memkdfields, cells; textcolor=:black)
overall = Plots.plot(a,b)
plots = plotsdir.("fields", "goalfield" .* ["cue", "mem", "cuemem"] .* "kde_vel>0")
sf(a, plots[1])
sf(b, plots[2])
sf(overall, plots[3])
a=field.show_all_fields(cuefields, cells; textcolor=:black)
b=field.show_all_fields(memfields, cells; textcolor=:black)
overall = Plots.plot(a,b)
plots = plotsdir.("fields", "goalfield" .* ["cue", "mem", "cuemem"] .* "hist_vel>0")
sf(a, plots[1])
sf(b, plots[2])
sf(overall, plots[3])

