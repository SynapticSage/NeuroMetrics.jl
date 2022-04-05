
## ,---.          |        ,---.o     |        |
# |  _.,---.,---.|        |__. .,---.|    ,---|
# |   ||   |,---||        |    ||---'|    |   |
# `---'`---'`---^`---'    `    ``---'`---'`---'
goalprops = ["currentPathLength","currentAngle"];
filters = merge(kws.filters,
                filt.correct,
                filt.notnan("currentAngle"), 
                filt.merge(filt.notnan("currentPathLength"), filt.minmax("currentPathLength", 2, 150)))
newkws = (;kws..., resolution=60, gaussian=3*0.5, filters=filters)
goalpath = field.get_fields(beh, spikes; props=goalprops, newkws...)
F["goal"] = goalpath
place_undergoal = field.get_fields(beh, spikes; props=["x","y"], newkws...)
F["place_undergoal"] = place_undergoal

if ploton
    a=field.plot.show_fields(goalpath.Rₕ)
    b=field.plot.show_fields(goalpath.Rₖ)
    grouping = field.group([goalpath.Rₕ,goalpath.Rₖ],["hist","kd"])
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

    P["goal"]            = model.run(spikes, beh, goalpath, props,
                                     type="goal")
    P["place_undergoal"] = model.run(spikes, beh, place_undergoal, props,
                                     type="place_undergoal")

    C["goal-place_undergoal"] = model.comparison.binaryop(P,
                                                          ["place_undergoal",
                                                           "goal"], op=-)
    C["goal/place_undergoal"] = model.comparison.binaryop(P,
                                                          ["place_undergoal",
                                                           "goal"], op=./)
    combine(groupby(C["goal-place_undergoal"], :area), :prob=>median)
    combine(groupby(C["goal/place_undergoal"], :area), :prob=>median)
    combine(groupby(C["goal-place_undergoal"], :area), :prob=>x->mean(x .> 0))
    combine(groupby(C["goal/place_undergoal"], :area), :prob=>x->mean(x .> 1))

    if ploton
        model.plot.individual_poisson_mean_unitarea(P["goal"])
        model.plot.individual_poisson_mean_unitarea(P["place_undergoal"])
        p = model.plot.comparison_binaryop_striplot_unitarea(C["goal-place_undergoal"],
                                                         hline=[0],
                                                         label="goal-place")
        p = model.plot.comparison_binaryop_striplot_unitarea(C["goal/place_undergoal"],
                                                         hline=[1],
                                                         label="goal/place")

        p = @df C["goal/place_undergoal"] density(:prob, group=:area)
        savefig(p, plotsdir("fields", "poisson_model", "goal-div-place-density.svg"))
        savefig(p, plotsdir("fields", "poisson_model", "goal-div-place-density.png"))
        savefig(p, plotsdir("fields", "poisson_model", "goal-div-place-density.pdf"))
        p = @df C["goal-place_undergoal"] density(:prob, group=:area)
        savefig(p, plotsdir("fields", "poisson_model", "goal-minus-place-density.svg"))
        savefig(p, plotsdir("fields", "poisson_model", "goal-minus-place-density.png"))
        savefig(p, plotsdir("fields", "poisson_model", "goal-minus-place-density.pdf"))
    end

end
