using Random, LinearAlgebra

abstract type OptimizationProblem end

function evaluate(::OptimizationProblem, genome::Vector{Float64}; kwargs...) :: Float64
    throw(MethodError(evaluate, (genome, kwargs...)))
end

function worse_than(::OptimizationProblem, first_fitness::Float64, second_fitness::Float64) :: Bool
    throw(MethodError(worse_than, (first_fitness, second_fitness)))
end

function problem_bounds(::OptimizationProblem) :: Bounds
    throw(MethodError(bounds, ()))
end

function maximize_fitness(::OptimizationProblem) :: Bool
    throw(MethodError(maximize_fitness, ()))
end

function equivalent(::OptimizationProblem, first_fitness::Float64, second_fitness::Float64) :: Bool
    return first_fitness == second_fitness
end

struct Bounds
    lower::Vector{Float64}
    upper::Vector{Float64}
end

"""
    FunctionProblem(; fitness_function, lower, upper, maximize=false)

Construct an optimization problem for a standard mathematical function.

# Arguments
- `fitness_function::Function`: The objective function to evaluate. It must accept a 
  single vector `x` and return a scalar numerical value.
- `lower::Vector`: A vector specifying the minimum allowable values for each dimension.
- `upper::Vector`: A vector specifying the maximum allowable values for each dimension.
- `maximize::Bool`: (Optional) Set to `true` for maximization, or `false` for 
  minimization (default).

# Examples
```julia
problem = FunctionProblem(
    fitness_function = x -> x[1]^2 + x[2]^2,
    lower = [-5.0, -5.0],
    upper = [5.0, 5.0]
)
```
"""
struct FunctionProblem <: OptimizationProblem
    fitness_function::Function
    _bounds::Bounds
    _maximize::Bool

    function FunctionProblem(
        fitness_function::Function,
        lower::Vector,
        upper::Vector,
        maximize::Bool
    )
        if length(lower) != length(upper)
            throw(ArgumentError("Lower and upper bound vectors must be the same length. Got $(length(lower)) and $(length(upper))."))
        end

        lower_f = Float64.(lower)
        upper_f = Float64.(upper)
        bounds = Bounds(lower_f, upper_f)
        new(fitness_function, bounds, maximize)
    end
end

function FunctionProblem(; fitness_function::Function, lower::Vector, upper::Vector, maximize::Bool=false)
    return FunctionProblem(fitness_function, lower, upper, maximize)
end

function evaluate(problem::FunctionProblem, genome::Vector{Float64}; kwargs...) :: Float64
   
    result = problem.fitness_function(genome; kwargs...)
    return result
end

function worse_than(problem::FunctionProblem, first_fitness::Float64, second_fitness::Float64)
    if isnan(first_fitness)
        return isnan(second_fitness) ? rand(Bool) : true
    elseif isnan(second_fitness)
        return false
    end
    return problem._maximize ? (first_fitness < second_fitness) : (first_fitness > second_fitness)
end

problem_bounds(problem::FunctionProblem) = problem._bounds
maximize_fitness(problem::FunctionProblem) = problem._maximize
