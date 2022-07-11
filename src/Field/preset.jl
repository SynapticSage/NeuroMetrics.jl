module preset

    import Field
    import Field: adaptive, fixed, metrics

    export field_presets, return_preset_funcs

    field_presets = Dict(

        :fixed=>
        (;fullname="fixed",
         desc=s"""classic fixed receptive field sampling""",
         fieldfunc   = fixed.RFs,
         gridfunc    = fixed.get_grid,
         occfunc     = fixed.get_occupancy,
         metricfuncs = metrics.information,
         postfunc    = nothing,
         grid_kws    = (;width=4)
        ),
        

        :legacy =>
        (;
        ),

        :yartsev => 
        (;fullname="adaptive-sampling-yartsev",
         desc=s"""adaptive sampling, yartsev 3d place field paper""",
         fieldfunc   = adaptive.yartsev,
         gridfunc    = adaptive.get_grid,
         occfunc     = adaptive.get_occupancy,
         metricfuncs = [metrics.information],
         postfunc    = nothing,
         grid_kws    = (;width=4)
        ),

        :skaggs => 
        (;fullname="adaptive sampling-skaggs",
         desc=s""" """,
         fieldfunc   = nothing,
         gridfunc    = nothing,
         occfunc     = nothing,
         metricfuncs = nothing,
         postfunc    = nothing,
        ),

       )


    function return_preset_funcs(p::Symbol)
        Tuple(field_presets[p][func] for func in [:fieldfunc, :gridfunc, 
                                                  :occfunc, :metricfuncs, 
                                                  :postfunc])
    end

end
