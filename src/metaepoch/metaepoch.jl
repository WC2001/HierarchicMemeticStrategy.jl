"""
    MetaepochRunner

An abstract type serving as the interface for all HMS search strategies.

# Interface
Any subtype of `MetaepochRunner` must implement the `run_metaepoch` function with the following signature:

```julia
run_metaepoch(
    runner::MetaepochRunner,
    fitness_function::Function,
    bounds::Bounds,
    initial_population::Vector{Vector{Float64}},
    minimize::Bool,
    options::Dict{String, Any}
)::MetaepochResult
```
"""
abstract type MetaepochRunner end

"""
    MetaepochResult(solution, best_fitness, populations, fitness_values, context)

A data structure that captures the outcome and trace of a single metaepoch execution.

# Fields
- `solution::Vector{Float64}`: The best input vector found during this metaepoch.
- `best_fitness::Float64`: The scalar fitness value associated with the `solution`.
- `populations::Vector{Vector{Vector{Float64}}}`: A history of populations. 
- `fitness_values::Vector{Vector{Float64}}`: A history of fitness scores. 

# Description
This object is returned by `run_metaepoch` and is used by the main HMS solver to 
update the tree state, check stop conditions, and build the final result.
"""
struct MetaepochResult
    solution::Vector{Float64}
    best_fitness::Float64
    populations::Vector{Vector{Vector{Float64}}}
    fitness_values::Vector{Vector{Float64}}
    context::Any
end

function run_metaepoch(
    ::MetaepochRunner,
    fitness_function::Function,
    bounds::Bounds,
    initial_population::Vector{Vector{Float64}},
    minimize::Bool,
    options::Dict{String,Any}
)::MetaepochResult
    error("run_metaepoch not implemented.")
end

