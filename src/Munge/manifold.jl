module manifold
    
    export get_dim_subset
    function get_dim_subset(dict_of_embeds::Dict, dim::Int)::Dict
        Dict(k=>v for (k,v) in dict_of_embeds
             if k[2] == dim)
    end
end

