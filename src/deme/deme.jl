using UUIDs
using Optim

mutable struct Deme
    id::String
    population::Population
    fitness_values::Union{Vector{Float64}, Nothing}
    level::Int
    best_fitness::Float64
    best_solution::Vector{Float64}
    best_solution_per_metaepoch::Vector{Any}
    best_fitness_per_metaepoch::Vector{Any}
    sprout::Union{Vector{Float64}, Nothing}
    parent_id::Union{String, Nothing}
    evaluations_count::Int
    is_active::Bool

    function Deme(
        problem::HMSOptimizationProblem,
        parent::Union{Deme, Nothing},
        population_size::Int,
        create_population::AbstractPopulationCreator,
        sigma::Vector{Vector{Float64}}
    )
        level = isnothing(parent) ? 1 : parent.level + 1
        sprout = isnothing(parent) ? get_u0(problem) : parent.best_solution_per_metaepoch[end]
        new_population_vector = create_population(;
            mean = sprout,
            lower = problem_bounds(problem).lower,
            upper = problem_bounds(problem).upper,
            population_size = population_size,
            tree_level = level,
            sigma = sigma
        )

        new_population = Population(new_population_vector, Float64[], problem)

        n_individuals = length(new_population.genomes)
        n_dimensions = length(new_population.genomes[1])
        if n_individuals != population_size || n_dimensions != length(problem_bounds(problem).lower)
            error("Created population is invalid.")
        end

        new(
            string(uuid4()),                         # id
            new_population,                          # population
            nothing,                                  # fitness_values
            level,                                    # level
            maximize_fitness(problem) ? -Inf : Inf,           # best_fitness (set to Inf for minimization)
            Float64[],                                # best_solution
            Any[],                                    # best_solution_per_metaepoch
            Any[],                                    # best_fitness_per_metaepoch
            sprout,                                   # sprout
            isnothing(parent) ? nothing : parent.id, # parent_id
            0,                                        # evaluations_count
            true                                      # is_active
        )
    end

    function Deme(;
        parent_id=nothing,
        level=1,
        best_solution=Float64[],
        best_solution_per_metaepoch=Any[], 
        best_fitness_values=Any[], 
        best_fitness=0.0, 
        evaluations_count=0,
        is_active=true
        )
        new(
            string(uuid4()),          # id
            Population(),             # population (empty)
            nothing,                  # fitness_values
            level,                         # level
            best_fitness,              # best_fitness
            best_solution,               # best_solution
            best_solution_per_metaepoch, # best_solution_per_metaepoch
            best_fitness_values,       # best_fitness_per_metaepoch
            nothing,                   # sprout
            parent_id,                 # parent_id
            evaluations_count,         # evaluations_count
            is_active                  # is_active
        )
    end

end

function update_deme!(deme::Deme, metaepoch_result::MetaepochResult, minimize::Bool = false)
    
    solution = metaepoch_result.solution
    value = metaepoch_result.best_fitness
    population = metaepoch_result.populations[end]

    deme.population.genomes = population
    push!(deme.best_fitness_per_metaepoch, value)
    push!(deme.best_solution_per_metaepoch, solution)

    if metaepoch_result.fitness_values !== nothing
        deme.fitness_values = metaepoch_result.fitness_values[end]
    end

    if metaepoch_result.context !== nothing
        deme.context = metaepoch_result.context
    end

    current_best = isempty(deme.best_fitness_per_metaepoch) ? (minimize ? Inf : -Inf) : deme.best_fitness
    is_better = minimize ? (value < current_best) : (value > current_best)
    
    if is_better
        deme.best_fitness = value
        deme.best_solution = solution
    end

    return deme
end

function update_deme!(deme::Deme, result::Optim.OptimizationResults, minimize::Bool = true)
    solution = Optim.minimizer(result)
    value = Optim.minimum(result)
    fval = minimize ? value : -value
    is_better = minimize ? (fval < deme.best_fitness) : (fval > deme.best_fitness)
    
    if is_better
        deme.best_fitness = fval
        deme.best_solution = solution
        add_individual(deme.population, solution, value)
    end

    return deme
end

function is_root(deme::Deme)::Bool
    return deme.parent_id === nothing
end

function get_not_processed_demes(demes::Vector{Deme}, processed_demes::Vector{Deme})::Vector{Deme}

    processed_ids = Set(deme.id for deme in processed_demes)
    return filter(deme -> !(deme.id in processed_ids), demes)
end

function get_leaves(demes::Vector{Deme}, max_height::Int)
    return filter(deme ->deme.level == max_height && length(deme.best_solution) > 0, demes)
end

struct FitnessResult
    fitness::Float64
    solution::Vector{Float64}
end

function get_best_solution(demes::Vector{Deme}, minimize::Bool)::FitnessResult
    op = minimize ? (<) : (>)
    best_fitness = minimize ? Inf : -Inf
    best_solution = nothing

    for deme in demes
        if isempty(deme.best_fitness)
            continue
        end

        if op(deme.best_fitness, best_fitness)
            best_fitness = deme.best_fitness
            best_solution = deme.best_solution
        end
    end

    return FitnessResult(best_fitness, best_solution)
end
