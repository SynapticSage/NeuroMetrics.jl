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
grouping = field.group([place.hist,place.kde],["hist","kd"])

if ploton
    name=plotsdir("fields", "individual_pf_vel=gt4")
    if !(isdir(name))
        mkdir(name)
    end
    name=plotsdir("fields", "individual_pf_vel=gt4", "place")
    individualgroups = field.plot.show_fieldgroups(grouping, as=Dict, nplots_factor=0.4)
    field.plot.save_dict_fieldplots(individualgroups, name)
end

if dopoisson

end
