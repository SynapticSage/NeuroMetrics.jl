
function remove_key_item(k::Dict, item)
    if item ∈ keys(k)
        pop!(k, item)
    end
    k
end

function lambda_keys(d::Dict, lambda::Function)
    d = copy(d)
    lambda_keys!(d, lambda)
end

function lambda_keys!(d::Dict, lambda::Function)
    for key ∈ keys(d)
        v = pop!(d, key)
        key = lambda(key)
        d[key] = v
    end
    return d
end
