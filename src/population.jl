using Statistics, Random, Distributions

mutable struct Population
    genomes::Vector{Vector{Float64}}
    fitnesses::Vector{Float64}
    problem::HMSOptimizationProblem

    function Population(genomes::Vector{Vector{Float64}}, fitnesses::Vector{Float64}, problem::HMSOptimizationProblem)
        new([Float64.(genome) for genome in genomes], Float64.(fitnesses), problem)
    end

    Population() = Population(Vector{Vector{Float64}}(), Float64[], FunctionProblem(fitness_function = x->x, lower = [0.0], upper = [1.0]))

end

function add_individual(pop::Population, genome::Vector{Float64}, fitness::Float64)
    push!(pop.genomes, genome)
    push!(pop.fitnesses, fitness)
end

function centroid(pop::Population)
    return Statistics.mean(reduce(hcat, pop.genomes), dims=2) |> vec
end


"""
    AbstractPopulationCreator

The abstract supertype for all population initialization strategies in HMS.

To implement a custom strategy, subtype this and implement the functor interface:
`(pc::MyCreator)(; mean, lower, upper, population_size, tree_level, sigma, kwargs...)::Vector{Vector{Float64}}`
"""

abstract type AbstractPopulationCreator end

"""
    reseed!(pc::AbstractPopulationCreator, rng::AbstractRNG)

Inject a new random number generator into the creator. 
This is used by the `hms` solver to ensure reproducibility when a `seed` is provided.
"""
reseed!(pc::AbstractPopulationCreator, rng::AbstractRNG) = pc

"""
    (pc::AbstractPopulationCreator)(; mean, lower, upper, population_size, tree_level, sigma, kwargs...)::Vector{Vector{Float64}}

Interface for creating populations. All subtypes must implement this keyword-based call signature.
"""
function (pc::AbstractPopulationCreator)(;
    mean::Union{Vector{Float64}, Nothing}, 
    lower::Vector{Float64}, 
    upper::Vector{Float64}, 
    population_size::Int, 
    tree_level::Int, 
    sigma::Vector{Vector{Float64}},
    kwargs...
)::Vector{Vector{Float64}}
    error("Subtype $(typeof(pc)) must implement the functional interface: (pc::$(typeof(pc)))(; mean, lower, upper, population_size, tree_level, sigma, kwargs...)::Vector{Vector{Float64}")
end


"""
    DefaultPopulationCreator()

The standard HMS initialization strategy.

### Behavior:
- **Level 1 (Root):** Uses a **Uniform distribution** across the entire search space. 
  This ensures global exploration at the top of the tree.
- **Level > 1 (Deme):** Uses a **Truncated Normal distribution** centered at the 
  `mean` (parent's best solution). The spread is controlled by `sigma[tree_level]`.


"""
mutable struct DefaultPopulationCreator <: AbstractPopulationCreator
    rng::AbstractRNG

    DefaultPopulationCreator() = new(Random.default_rng())
end


function reseed!(pc::DefaultPopulationCreator, rng::AbstractRNG)
    pc.rng = rng
    return pc
end

function (pc::DefaultPopulationCreator)(;
    mean::Union{Vector{Float64}, Nothing}, 
    lower::Vector{Float64}, 
    upper::Vector{Float64}, 
    population_size::Int, 
    tree_level::Int, 
    sigma::Vector{Vector{Float64}},
    kwargs...
)::Vector{Vector{Float64}}
    if tree_level == 1
        return [
            [rand(pc.rng, Uniform(lower[i], upper[i])) for i in eachindex(lower)] 
            for _ in 1:population_size
        ]
    else
        sd = sigma[tree_level]
        genomes = Vector{Vector{Float64}}()
        
        for _ in 1:(population_size - 1)
            individual = [
                rand(pc.rng, Truncated(Normal(mean[i], sd[i]), lower[i], upper[i]))
                for i in eachindex(mean)
            ]
            push!(genomes, individual)
        end
        
        push!(genomes, copy(mean))
        return genomes
    end
end

"""
    NormalPopulationCreator()

A focused initialization strategy that uses Normal distributions for all levels.

### Behavior:
- Uses a **Truncated Normal distribution** for every tree level.
- If no `mean` is provided at Level 1, it centers the population at the search space midpoint.
- Ideal for optimization where a good initial guess is available.
"""
mutable struct NormalPopulationCreator <: AbstractPopulationCreator
    rng::AbstractRNG
    
    NormalPopulationCreator() = new(Random.default_rng())
end

function reseed!(pc::NormalPopulationCreator, rng::AbstractRNG)
    pc.rng = rng
    return pc
end

function (pc::NormalPopulationCreator)(;
    mean::Union{Vector{Float64}, Nothing}, 
    lower::Vector{Float64}, 
    upper::Vector{Float64}, 
    population_size::Int, 
    tree_level::Int, 
    sigma::Vector{Vector{Float64}},
    kwargs...
)::Vector{Vector{Float64}}

    mean = if tree_level == 1 && isnothing(mean)
        (lower .+ upper) ./ 2.0
    else
        mean
    end

    sd = sigma[tree_level]
    genomes = Vector{Vector{Float64}}()
    
    for _ in 1:(population_size - 1)
        individual = [
            rand(pc.rng, Truncated(Normal(mean[i], sd[i]), lower[i], upper[i]))
            for i in eachindex(mean)
        ]
        push!(genomes, individual)
    end
    
    push!(genomes, copy(mean))
    return genomes
end