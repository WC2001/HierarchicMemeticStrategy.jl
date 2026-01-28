euclidean(a::Vector{Float64}, b::Vector{Float64}) = sqrt(sum((a .- b).^2))
manhattan(a::Vector{Float64}, b::Vector{Float64}) = sum(abs.(a .- b))

function sprout_default_euclidean_distances(sigma::Vector{Vector{Float64}})
    sprouting_condition_distance_ratio = 0.6
    return [sum(s .* sprouting_condition_distance_ratio) for s in sigma]
end

"""
    AbstractSproutCondition

The abstract supertype for all sprout validation strategies in HMS.

Sprout conditions determine whether a new deme should be created at a specific 
location. This is typically used to maintain diversity by preventing demes from 
clustering too closely together.

### To Implement a Custom Strategy:
Subtype this and implement the keyword-based functor interface:
`(sc::MyCondition)(potential_sprout, level, demes; kwargs...)::Bool`
"""
abstract type AbstractSproutCondition end

"""
    (sc::AbstractSproutCondition)(potential_sprout, potential_sprout_level, demes; kwargs...)::Bool

The core logic for sprout condition. Returns `true` if a new deme is allowed to sprout.
- `potential_sprout`: The sprout to create new deme from.
- `potential_sprout_level`: The level in HMS tree where the deme would be created.
- `demes`: A list of all currently active demes in the HMS tree.
- `kwargs`: Additional context.
"""
function (sc::AbstractSproutCondition)(
    potential_sprout::Vector{Float64}, 
    potential_sprout_level::Int, 
    demes::Vector{Deme}; 
    kwargs...
)::Bool
    error("Subtype $(typeof(sc)) must implement the functor interface.")
end


"""
    MetricSproutCondition(metric::Function, max_distances::Vector{Float64})

A distance-based validation strategy that uses a metric to ensure demes 
stay separated.

### Fields
- `metric`: A function `(a, b) -> Float64` (e.g., `euclidean` or `manhattan`).
- `max_distances`: A vector of thresholds. `max_distances[i]` is the minimum 
  distance required between a new sprout at level `i+1` and existing demes at that level.

### Behavior
A sprout is rejected if the distance to the centroid of **any** existing deme 
at the target level is **less than or equal to** the specified threshold.
"""
struct MetricSproutCondition <: AbstractSproutCondition
    metric::Function
    max_distances::Vector{Float64}
end

function (sc::MetricSproutCondition)(
    potential_sprout::Vector{Float64}, 
    potential_sprout_level::Int, 
    demes::Vector{Deme}; 
    kwargs...
)::Bool
    level_demes = filter(d -> d.level == potential_sprout_level, demes)
    
    threshold = sc.max_distances[potential_sprout_level - 1]

    for deme in level_demes
        if sc.metric(centroid(deme.population), potential_sprout) <= threshold
            return false 
        end
    end
    
    return true
end


"""
    DefaultSproutCondition(sigma::Vector{Vector{Float64}})

The default HMS sprout condition using Euclidean distance.
"""
function DefaultSproutCondition(sigma::Vector{Vector{Float64}})
    return MetricSproutCondition(euclidean, sprout_default_euclidean_distances(sigma))
end
