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


