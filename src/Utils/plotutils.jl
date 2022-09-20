module plotutils

    import ..Utils
    import Plots
    using Dates
    using ColorSchemes

    export  plotcolor


    function set_theme_timebased(time::Union{Real,Nothing}=nothing)
        time = time === nothing ? hour(now()) : time
        if Utils.in_range(time, [0,5]) ||
           Utils.in_range(time, [20, 24])
            Plots.theme(:brigh)
            theme="dark"
        else
            Plots.theme(:bright)
            theme="bright"
        end
    end
    
    function plotcolor(propvec::AbstractVector, cmap::Symbol)
        propvec = Utils.norm_extrema(propvec)
        get.([colorschemes[cmap]], propvec)
    end

end