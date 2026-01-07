# Abstract type for genome types
abstract type GenomeType end

# Concrete genome types
struct RealGenome <: GenomeType end
struct BinaryGenome <: GenomeType end
struct PermutationGenome <: GenomeType end

# Individual struct with parameterized genome type
mutable struct Individual{T<:GenomeType}
    genome::Vector{Float64}
    fitness::Float64
    problem::OptimizationProblem
    genome_type::T

    function Individual{T}(genome::Vector{<:Real}, problem::OptimizationProblem, genome_type::T, fitness::Real = NaN,) where {T<:GenomeType}
        new{T}(Float64.(genome), Float64(fitness), problem, genome_type)
    end
end
