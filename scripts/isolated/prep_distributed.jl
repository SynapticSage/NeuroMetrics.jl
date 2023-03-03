using Distributed
using DrWatson

while length(workers()) < 16
    try
        addprocs(2)
        println("Workers: $(length(workers))")
    catch exception
        print(exception)
    end
end

@everywhere begin 
    using DrWatson
    quickactivate(expanduser("~/Projects/goal-code"))
    include(scriptsdir("isolated", "imports_isolated.jl"))
end

include(scriptsdir("isolated","load_cyclewise_checkpoint.jl"))
