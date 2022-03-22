module utils

import Random

# FIND COLOR SCHEMES (using matlab utilties)
#settings = begin
#    function set_name(settings, name; N=10, method="cmocean")
#        ref = settings["colornames"][name]
#        if method == "cmocean"
#            mat"$ref = cmocean($ref, $N)";
#        elseif method == "crameri"
#            mat"$ref = crameri($ref, $N)";
#        end
#        settings["colordict"][name] = ref;
#        return settings
#    end
#    settings = set_name(settings, "hpc");
#    settings = set_name(settings, "pfc");
#    settings = set_name(settings, "all", N=maximum(downsamp.unit));
#    settings = set_name(settings, "all_crameri",
#                        method="crameri",
#                        N=maximum(downsamp.unit))
#end

    #                            
    # .   .|    o|    o|         
    # |   ||--- .|    .|--- ,   .
    # |   ||    ||    ||    |   |
    # `---'`---'``---'``---'`---|
    #                       `---'
    function randomize_int(X)
        Xmin = minimum(X);
        Xmax = maximum(X);
        initial = collect(Xmin:Xmax);
        final   = Random.shuffle(initial);
        mapping(x) = Dict(initial .=> final)[x]
        map(mapping, X)
    end
end
