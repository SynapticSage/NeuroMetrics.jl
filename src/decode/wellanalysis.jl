using StatsBase
export cosine_similarity_to_well

"""
cosine_similarity_to_well
gets cosine similarity of a given decode's vector to the vector between
the animal and a well. This is to observe the evolution of sequence vectors
between an animal's decodes and his wells
"""
function cosine_similarity_to_well(X, well; decode_vec_method=:decode,
        unit_decode=false)
    unitvec(x⃗) = x⃗ ./ abs.(x⃗)
    vector_animaltowell = (well.x .- X.start_x) .- (well.y .- X.start_y)im;
    vector_animaltowell = unitvec(vector_animaltowell);
    if decode_vec_method == :decode
        vector_decode = (X.stop_x_dec .- X.start_x_dec) .+ (X.stop_y_dec .- X.start_y_dec)im
    elseif decode_vec_method == :animal_to_decode_end
        vector_decode = (X.stop_x_dec .- X.start_x)     .+ (X.stop_y_dec .- X.start_y)im
    end
    if unit_decode
        vector_decode = unitvec(vector_decode);
        θ = angle.(vector_decode .- vector_animaltowell)
    else
        θ = angle.(unitvec(vector_decode) .- vector_animaltowell)
    end
    abs.(vector_decode) .* abs.(vector_animaltowell) .* cos.(θ);
end

# Cosine similarity to future₁, future₂, past₁, past₂ 
# Split for (home, arena), (correct, incorrect)
# --------------------------------------------------
# (Accomplished by getting a table for each of these 4, and averaging a
# categorical symbolizing the two split tuples)

function cosine_similarity_metric(X::DataFrame, beh::DataFrame,
        well::DataFrame; kws...)
    homewell   = fit(Histogram, beh.stopWell)
    beh.ishome = beh.stopWell == homewell
    splits = [:ishome, :correct]
    register(X, beh, )
    groupby(beh, splits)
end
