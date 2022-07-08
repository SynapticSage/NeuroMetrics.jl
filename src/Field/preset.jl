module preset

    import Field
    import Field: adaptive


    field_presets = Dict(

        :fixed,

        :legacy,

        :adaptive => 
        (;fullname="",
         desc=s""" """,
         fieldfunc=nothing,
         postfunc=nothing,
         gridfunc=nothing,
         occfunc=nothing,
         metricfuncs=nothing,
        ),
       )
end
