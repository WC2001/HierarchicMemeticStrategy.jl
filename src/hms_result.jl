"""
    HMSResult

A structure containing the results of an HMS optimization run.
"""
struct HMSResult
    metaepoch_summaries::Vector{MetaepochSummary}
    hms_result_visualizer::HMSResultVisualizer
    optimization_problem::HMSOptimizationProblem
    solution::Vector{Float64}
    best_fitness::Float64
    fitness_evaluation_count::Int
end

function HMSResult(metaepoch_summaries::Vector{MetaepochSummary}, hms_result_visualizer::HMSResultVisualizer, problem::HMSOptimizationProblem)
    last_summary = metaepoch_summaries[end]
    solution = last_summary.solution
    best_fitness = last_summary.best_fitness
    fitness_evaluation_count = last_summary.fitness_evaluation_count

    HMSResult(metaepoch_summaries, hms_result_visualizer, problem, solution, best_fitness, fitness_evaluation_count)
end

function Base.show(io::IO, result::HMSResult)
    println(io, "┌────────────────────── HMSResult Summary ───────────────────────")
    println(io, "│ • Metaepoch Count:      ", length(result.metaepoch_summaries))
    println(io, "│ • Solution:             ", result.solution)
    println(io, "│ • Best Fitness:         ", result.best_fitness)
    println(io, "│ • Fitness Evaluations:  ", result.fitness_evaluation_count)
    println(io, "└────────────────────────────────────────────────────────────────")
end

"""
    summary(result::HMSResult)

Print a formatted summary table to the standard output for a quick overview of the results.
"""
function summary(result::HMSResult)
    show(stdout, result)
end

"""
    f_calls(result::HMSResult)

Return the total number of times the objective function was evaluated across all demes.
"""
function f_calls(result::HMSResult)
    return result.fitness_evaluation_count
end

"""
    iterations(result::HMSResult)

Return the total count of metaepochs performed during the run.
"""
function iterations(result::HMSResult)
    return length(result.metaepoch_summaries)
end

"""
    solution(result::HMSResult)

Return the best-found input vector (solution).
"""
function solution(result::HMSResult)
    return result.solution
end

"""
    best_fitness(result::HMSResult)

Returns the scalar fitness value at the best-found point.
"""
function best_fitness(result::HMSResult)
    return result.best_fitness
end

"""
    metaepoch_data(result::HMSResult)

Return algorithm data.
"""
function metaepoch_data(result::HMSResult)
    return result.metaepoch_summaries
end

"""
    plotPopulations(result::HMSResult; x_index=1, y_index=2, filename="")

Show interactive plot presenting changes in populations during metaepochs.

# Arguments
- `result::HMSResult`: The result object returned by the optimizer.
- `x_index::Int`: The index of the dimension to plot on the x-axis (default: `1`).
- `y_index::Int`: The index of the dimension to plot on the y-axis (default: `2`).
- `filename::String`: (Optional) The name of the output HTML file. If provided without the `.html` 
  extension, it will be added automatically. If left empty, a default name with a timestamp 
  is generated (e.g., `hms_plot_2026-01-11_12-00-00.html`).
"""
function plotPopulations(result::HMSResult; x_index::Int=1, y_index::Int=2, filename::String="")
    populations_over_metaepochs = plotPopulations(
        result.hms_result_visualizer,
        result.optimization_problem,
        x_index, 
        y_index,
        filename
    )
    display(populations_over_metaepochs)
end

"""
    plotBestFitness(result::HMSResult)

Plot best_fitness value during metaepochs.
"""
function plotBestFitness(result::HMSResult)
    best_fitness_over_metaepochs = plotBestFitness(result.hms_result_visualizer)
    display(best_fitness_over_metaepochs)
end

"""
     plotDeme(result::HMSResult; deme_index::Int=1, x_index::Int=1, y_index::Int=2)

Plot single Deme population over metaepochs.

# Arguments
- `result::HMSResult`: The result object returned by the HMS optimizer.
- `deme_index::Int`: The index of deme to visualize (default: `1`).
- `x_index::Int`: The index of the dimension to plot on the x-axis (default: `1`).
- `y_index::Int`: The index of the dimension to plot on the y-axis (default: `2`).
"""
function plotDeme(result::HMSResult; deme_index::Int=1, x_index::Int=1, y_index::Int=2)
    deme_plot = plotDemeHistory(result.hms_result_visualizer, result.optimization_problem, deme_index, x_index, y_index)
    display(deme_plot)
end

