if :ploton ∉ propertynames(Main)
    ploton = dopoissonmodel = false
end

#                                                   
# ,---.|                       ,---.o     |        |
# |---'|    ,---.,---.,---.    |__. .,---.|    ,---|
# |    |    ,---||    |---'    |    ||---'|    |   |
# `    `---'`---^`---'`---'    `    ``---'`---'`---'
props = ["x", "y"]
newkws = (; kws..., props=props, filters=merge(kws.filters))
place = field.get_fields(beh, spikes; newkws...);
F["place"] = place

# Generate all of the fields >  5cm/s
# --------------------------
if ploton
    grouping = field.group([place.Rₕ,place.Rₖ],["hist","kd"])
    a        = field.plot.show_fields(place.Rₕ)
    b        = field.plot.show_fields(place.Rₖ)
    overall = field.plot.show_fieldgroups(grouping)
    name=plotsdir("fields", "individual_pf_vel=gt4")
    if !(isdir(name))
        mkdir(name)
    end
    name=plotsdir("fields", "individual_pf_vel=gt4", "place")
    individualgroups = field.plot.show_fieldgroups(grouping, as=Dict,
                                                   nplots_factor=0.4)
    field.plot.save_dict_fieldplots(individualgroups, name)
end

if dopoissonmodel
    # Ready the data
    data = field.model.data(spikes, beh; grid=place.cgrid, props=props)
    # Acquire the probabilities fields | data
    likelihood = field.model.probability(data, place.Rₕ)
    # Convert to dataframe
    likelihood = table.to_dataframe(likelihood, other_labels=Dict(:type=>"xy"),
                                    name="prob")
    likelihood[!,"logprob"] = log10.(likelihood.prob)
    table.naninf_to_missing!(likelihood, [:prob, :logprob])
    P["place"] = likelihood

    if ploton
        model.plot.individual_poisson_mean_unitarea(likelihood)
    end
end
