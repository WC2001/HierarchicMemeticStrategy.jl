using Random, LinearAlgebra

abstract type HMSOptimizationProblem end

function get_u0(problem::HMSOptimizationProblem)
    return hasfield(typeof(problem), :u0) ? problem.u0 : nothing
end

function evaluate(::HMSOptimizationProblem, genome::Vector{Float64}; kwargs...) :: Float64
    throw(MethodError(evaluate, (genome, kwargs...)))
end

function worse_than(::HMSOptimizationProblem, first_fitness::Float64, second_fitness::Float64) :: Bool
    throw(MethodError(worse_than, (first_fitness, second_fitness)))
end

function problem_bounds(::HMSOptimizationProblem) :: Bounds
    throw(MethodError(bounds, ()))
end

function maximize_fitness(::HMSOptimizationProblem) :: Bool
    throw(MethodError(maximize_fitness, ()))
end

function equivalent(::HMSOptimizationProblem, first_fitness::Float64, second_fitness::Float64) :: Bool
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
- `u0::Union{Vector{Float64}, Nothing}`: (Optional) An initial guess or starting point. 
When provided, certain population creators (like `NormalPopulationCreator`) can use 
this as the center of the initial root population.

# Examples
```julia
problem = FunctionProblem(
    fitness_function = x -> x[1]^2 + x[2]^2,
    lower = [-5.0, -5.0],
    upper = [5.0, 5.0]
)
```
"""
struct FunctionProblem <: HMSOptimizationProblem
    fitness_function::Function
    u0::Union{Vector{Float64}, Nothing}
    _bounds::Bounds
    _maximize::Bool

    function FunctionProblem(
        fitness_function::Function,
        lower::Vector,
        upper::Vector,
        maximize::Bool,
        u0::Union{Vector{Float64}, Nothing}
    )
        if length(lower) != length(upper)
            throw(ArgumentError("Lower and upper bound vectors must be the same length. Got $(length(lower)) and $(length(upper))."))
        end

        lower_f = Float64.(lower)
        upper_f = Float64.(upper)
        bounds = Bounds(lower_f, upper_f)
        new(fitness_function, u0, bounds, maximize)
    end
end

function FunctionProblem(;
    fitness_function::Function, 
    lower::Vector,
    upper::Vector, 
    maximize::Bool=false, 
    u0::Union{Vector{Float64}, Nothing}=nothing
)
    return FunctionProblem(fitness_function, lower, upper, maximize, u0)
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
