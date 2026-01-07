euclidean(a::Vector{Float64}, b::Vector{Float64}) = sqrt(sum((a .- b).^2))

function sprout_default_euclidean_distances(sigma::Vector{Vector{Float64}})
    sprouting_condition_distance_ratio = 0.6
    return [sum(s .* sprouting_condition_distance_ratio) for s in sigma]
end


"""
    sc_max_metric(metric::Function, max_distances::Vector{Float64})

Create a sprouting condition that ensures a new deme is sufficiently far from existing demes at the same level.

# Arguments
- `metric::Function`: A distance function, e.g., `dist(a, b)`, used to calculate the separation between points.
- `max_distances::Vector{Float64}`: A vector where `max_distances[i]` defines the minimum required distance for level `i + 1`.

# Returns
- A function `(potential_sprout, potential_sprout_level, demes) -> Bool` that returns `true` if the potential sprout location is valid.

# Description
If the distance to **any** existing deme at that level is less than the threshold defined in `max_distances`, the condition returns `false`, and the sprout is suppressed. This mechanism is critical for:
1. **Diversity**: Preventing multiple demes from converging on the same local optimum.
2. **Exploration**: Forcing the algorithm to seek out new, unexplored regions of the search space.

# Example
```julia
# Use Euclidean distance with a 2.0 threshold for level 2 and 0.5 for level 3
euclidean(a::Vector{Float64}, b::Vector{Float64}) = sqrt(sum((a .- b).^2))
sprout_cond = sc_max_metric(euclidean, [2.0, 0.5])
```
"""
function sc_max_metric(metric::Function, max_distances::Vector{Float64})
    return function(potential_sprout::Vector{Float64}, potential_sprout_level::Int, demes::Vector{Deme})
        level_demes = filter(d -> d.level == potential_sprout_level, demes)

        function single_deme_condition(deme)
            deme_centroid = centroid(deme.population)
            return metric(deme_centroid, potential_sprout) > max_distances[potential_sprout_level - 1]
        end

        return all(deme -> single_deme_condition(deme), level_demes)
    end
end

const DEFAULT_SC = function (sigma::Vector{Vector{Float64}})
    return sc_max_metric(euclidean, sprout_default_euclidean_distances(sigma))
    
end


manhattan(a::Vector{Float64}, b::Vector{Float64}) = sum(abs.(a .- b))