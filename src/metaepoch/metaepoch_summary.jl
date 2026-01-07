mutable struct MetaepochSummary
    demes::Vector{Deme}
    best_fitness::Float64
    solution::Vector{Float64}
    time::Float64
    fitness_evaluation_count::Int
    blocked_sprouts::Vector{Vector{Float64}}
end


MetaepochSummary(; fitness_evaluation_count::Int=0, demes::Vector{Deme}=Deme[]) = MetaepochSummary(
    demes,
    0.0,
    Float64[],
    0.0,
    fitness_evaluation_count,
    Vector{Vector{Float64}}()
)

function log_metaepoch_summary(
    summary::MetaepochSummary,
    metaepoch_count::Int,
    root::Deme,
    best::FitnessResult,
    log_tree::Bool
)
    if log_tree 
        visualizer = HMSTreeVisualizer(root, summary.demes, best)
        printTree(visualizer)
        println()
    end

    println("Metaepoch: ", metaepoch_count, "  Best fitness: ", summary.best_fitness)
    println("Best solution: ", summary.solution)
    valid_demes = count(d -> !isempty(d.best_solution), summary.demes)
    println("Demes count: ", valid_demes)
    println("Fitness evaluations: ", summary.fitness_evaluation_count)
    println("Blocked sprouts: ", summary.blocked_sprouts)
    println()

end
