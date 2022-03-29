
## ,---.          |        ,---.o     |        |
# |  _.,---.,---.|        |__. .,---.|    ,---|
# |   ||   |,---||        |    ||---'|    |   |
# `---'`---'`---^`---'    `    ``---'`---'`---'
goalprops = ["currentPathLength","currentAngle"];
filters = merge(kws.filters,
                filt.notnan("currentPathLength"), 
                filt.correct,
                filt.notnan("currentAngle"), 
                filt.minmax("currentPathLength", 2, 150))
newkws = (;kws..., resolution=80, gaussian=3*0.5, filters=filters)
goalpath = field.get_fields(beh, spikes; props=goalprops, newkws...)
F["goal"] = goalpath
place_undergoal = field.get_fields(beh, spikes; props=["x","y"], newkws...)
F["place_undergoal"] = place_undergoal

if ploton
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
end

sp = raw.filterTables(spikes, filters=filters)[1]
summary = combine(groupby(sp,:unit), nrow)
sum(summary.nrow .< 100)

if dopoissonmodel

    P["goal"]            = model.run(spikes, beh, goalpath,        goalprops)
    P["place_undergoal"] = model.run(spikes, beh, place_undergoal, props)

    if ploton
        model.plot.individual_poisson_mean_unitarea(P["goal"])
    end

end
