using Statistics

mutable struct Population
    genomes::Vector{Vector{Float64}}
    fitnesses::Vector{Float64}
    problem::OptimizationProblem

    function Population(genomes::Vector{Vector{Float64}}, fitnesses::Vector{Float64}, problem::OptimizationProblem)
        new([Float64.(genome) for genome in genomes], Float64.(fitnesses), problem)
    end

    Population() = Population(Vector{Vector{Float64}}(), Float64[], FunctionProblem(x->x,[0.0], [1.0], false))

end

function add_individual(pop::Population, genome::Vector{Float64}, fitness::Float64)
    push!(pop.genomes, genome)
    push!(pop.fitnesses, fitness)
end

function centroid(pop::Population)
    return Statistics.mean(reduce(hcat, pop.genomes), dims=2) |> vec
end