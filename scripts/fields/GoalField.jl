
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

    # Ready the data
    data = field.model.data(spikes, beh; grid=goalpath.cgrid, props=props)
    # Acquire the probabilities fields | data
    likelihood = field.model.probability(data, place.hist)
    # Convert to dataframe
    likelihood = field.to_dataframe(likelihood, other_labels=Dict(:type=>"goal"),
                                    name="prob")
    likelihood[!,"logprob"] = log10.(likelihood.prob)
    table.naninf_to_missing!(likelihood, [:prob, :logprob])
    P["goal"] = likelihood

    if ploton
        model.plot.individual_poisson_mean_unitarea(likelihood)
    end

end
