module GoalFetchAnalysis

    push!(LOAD_PATH, @__DIR__)

    __revise_mode__ = :eval
    using OhMyREPL

    # General
    using Load   # Loading/saving/manipulating raw data
    @info "got here"
    using Utils # General utilites
    @info "got here"
    using Utils # General utilites
    import Filt
    @info "got here"
    using Utils # General utilites
    @info "got here"
    import Shuffle # General receptive field codes
    @info "got here"

    import Table # Manipulating tables
    @info "got here"

    # Field related
    import Field # General receptive field codes
    @info "got here"
    import Timeshift
    @info "got here"

    # Goal vector measures
    import Statistic

    # Deode Related
    import Decode # General receptive field codes

    export Load
    export Utils
    export Table 
    export Field, Timeshift
    export Decode

end
