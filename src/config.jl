"""
    TreeLevelConfig(metaepoch_type, options)

Configuration for a specific level in the HMS hierarchy.

# Fields
- `metaepoch_type::Type{<:MetaepochRunner}`: The evolutionary strategy or algorithm 
  to run at this level (e.g., `EvolutionaryGAMetaepoch` or `EvolutionaryCMAESMetaepoch`).
- `options::Dict{String, Any}`: A dictionary of hyperparameters specific to the chosen 
  `metaepoch_type` (e.g. crossover rates, mutation).

# Usage
Each level of the HMS tree can employ a different search strategy. By providing a 
`Vector{TreeLevelConfig}` to the solver, you define how the search intensifies 
from the root level down to the leaves.

# Example
```julia
level_1 = TreeLevelConfig(EvolutionaryGAMetaepoch, Dict("crossoverRate" => 0.3, "mutationRate" => 0.5))
```
"""
mutable struct TreeLevelConfig
    metaepoch_type::Type{<:MetaepochRunner}   
    options::Dict{String, Any}

end

function TreeLevelConfig(metaepoch_type::Type{<:MetaepochRunner})
    TreeLevelConfig(metaepoch_type, Dict{String, Any}())
end

function set_option!(config::TreeLevelConfig, key::String, value::Any)
    config.options[key] = value
end

""" DEFAULT_LEVEL_CONFIG

The default HMS tree structure used by the HMS solver if no custom configuration is provided.

It defines a 2-level tree:

Level 1: Uses EvolutionaryGAMetaepoch for broad exploration.

Level 2: Uses EvolutionaryCMAESMetaepoch for refined local search. 
""" 
const DEFAULT_LEVEL_CONFIG::Vector{TreeLevelConfig} = [
    TreeLevelConfig(EvolutionaryGAMetaepoch),
    TreeLevelConfig(EvolutionaryCMAESMetaepoch)
]

function run_metaepoch_at_level(
    metaepoch_config::TreeLevelConfig,
    fitness_function::Function,
    bounds::Bounds,
    initial_population::Vector{Vector{Float64}},
    minimize::Bool,
    sigma::Vector{Float64} = nothing
)::MetaepochResult
    metaepoch_type = metaepoch_config.metaepoch_type
    options = metaepoch_config.options
    options["sigma"] = sigma
    if metaepoch_type <: MetaepochRunner
        return run_metaepoch(
            metaepoch_type(), 
            fitness_function, 
            bounds, 
            initial_population, 
            minimize, 
            options
        )
    else
        error("Unsupported metaepoch type: $metaepoch_type")
    end

end


"""
    default_population_sizes(tree_height::Int)

Generate set of population sizes for each level of the HMS tree.

# Arguments
- `tree_height::Int`: Total number of levels in the HMS hierarchy.

# Returns
- `Vector{Int}`: A vector of population sizes, one for each level, decreasing as the tree gets deeper.

# Description
In a Hierarchical Multi-Strategy search, it is often beneficial to have larger populations 
at the root level to maintain diversity and smaller populations at lower levels for 
specialized local search. This function implements a geometric decay:
- **Root Level**: Starts with a population of 60.
- **Decay**: Each subsequent level's population is reduced by 50%.

# Example
```julia
sizes = default_population_sizes(4)
# Returns: [60, 30, 15, 8]
```
"""
function default_population_sizes(tree_height::Int)
    population_size = 60
    population_size_ratio = 0.5
    population_sizes = Int[]

    for height in 1:tree_height
        push!(population_sizes, population_size)
        population_size = round(Int, population_size * population_size_ratio)
    end

    return population_sizes
end


"""
    default_sigma(lower::Vector{Float64}, upper::Vector{Float64}, tree_height::Int)

Generate sigma parameter for each level of the HMS tree.

# Arguments
- `lower::Vector{Float64}`: Lower bounds of the search space.
- `upper::Vector{Float64}`: Upper bounds of the search space.
- `tree_height::Int`: Total number of levels in the HMS hierarchy.

# Returns
- `Vector{Vector{Float64}}`: A list of sigma vectors, one for each level, where values decrease geometrically as the tree gets deeper.

# Description
Sigma value is a fraction of the total domain length (`upper - lower`). This function 
applies a geometric decay:
- **Root Level**: Starts at 4% of the domain length (`sigma_ratio = 0.04`).
- **Subsequent Levels**: Each level's sigma is reduced by 50% (`sigma_exponent = 0.5`) compared to the level above.

This ensures that top-level demes perform broad global exploration, while lower-level demes focus on increasingly fine-grained local exploitation.
"""
function default_sigma(lower::Vector{Float64}, upper::Vector{Float64}, tree_height::Int)
    sigma_ratio = 0.04
    sigma_exponent = 0.5
    domain_length = upper .- lower
    sigma = Vector{Vector{Float64}}()

    for height in 1:tree_height
        push!(sigma, domain_length .* sigma_ratio)
        sigma_ratio *= sigma_exponent
    end

    return sigma
end