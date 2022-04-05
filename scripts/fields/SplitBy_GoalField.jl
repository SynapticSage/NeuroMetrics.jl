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
if ploton
    name=plotsdir("fields", "individual_gfSplitBy_vel=gt2")
    mkdir(name)
    name=plotsdir("fields", "individual_gfSplitBy_vel=gt2", "goal")
    individualgroups = field.plot.show_fieldgroups(grouping, as=Dict,
                                                   nplots_factor=0.4);
    field.plot.save_dict_fieldplots(individualgroups, name,
                                    ext=["svg","png","pdf"])
end
