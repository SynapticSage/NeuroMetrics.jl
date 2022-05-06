module behavior
    function add_next_target(beh)
    end
    function add_previous_target(beh)
    end
    additions = Dict("previous_target" => add_previous_target,
                     "next_target"=> add_next_target)
end
export behavior

