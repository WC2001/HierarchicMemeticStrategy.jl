abstract type LocalStopCondition end

struct FitnessStable <: LocalStopCondition
    max_deviation::Float64
    n_metaepochs::Int
end

function (lsc::FitnessStable)(deme::Deme)::Bool
    if lsc.n_metaepochs > length(deme.best_fitness_per_metaepoch)
        return false
    end

    avg_fits = [
        mean([ind.fitness for generation in deme._history[n] for ind in generation])
        for n in -lsc.n_metaepochs:-1
    ]

    return mean(avg_fits) - minimum(avg_fits) <= lsc.max_deviation
end

"""
    LocalProblemEvaluationLimitReached(limit::Int)

A `LocalStopCondition` that terminates a specific deme once it has performed a 
certain number of objective function evaluations.

# Fields
- `limit::Int`: The maximum number of evaluations allowed for the specific deme.

# Example
```julia
# Limit each leaf deme to 1000 evaluations
lsc = LocalProblemEvaluationLimitReached(1000)
```
"""
struct LocalProblemEvaluationLimitReached <: LocalStopCondition
    limit::Int

    function LocalProblemEvaluationLimitReached(limit::Int)
        limit > 0 || throw(ArgumentError("limit must be a positive integer, got $limit"))
        new(limit)
    end
end

function(lsc::LocalProblemEvaluationLimitReached)(deme::Deme, metaepoch_summaries::Vector{MetaepochSummary})::Bool
    return deme.evaluations_count >= lsc.limit
end

"""
    AllChildrenStopped()

A `LocalStopCondition` that triggers only when all children of the current deme 
have reached their own stopping criteria.
"""
struct AllChildrenStopped <: LocalStopCondition end


function (lsc::AllChildrenStopped)(deme::Deme, metaepoch_summaries::Vector{MetaepochSummary})::Bool

    isempty(metaepoch_summaries) && return false

    last_metaepoch_summary = metaepoch_summaries[end]
    demes = last_metaepoch_summary.demes
    children = filter(d -> d.parent_id == deme.id, demes)

    isempty(children) && return false

    return all(!child.is_active for child in children)
end

"""
    MetaepochWithoutBestFitnessImprovement(n_metaepochs::Int)

A `LocalStopCondition` that triggers if the best fitness within a deme fails to 
improve for a specified number of consecutive metaepochs.

# Fields
- `n_metaepochs::Int`: The threshold of consecutive stagnant metaepochs allowed before stopping.

# Description
This condition is used to prune demes that have reached a local optimum or are 
stuck in a plateau.

# Example
```julia
# Stop a deme if it doesn't improve for 5 consecutive metaepochs
lsc = MetaepochWithoutBestFitnessImprovement(5)
```
"""
struct MetaepochWithoutBestFitnessImprovement <: LocalStopCondition 
    n_metaepochs::Int

    function MetaepochWithoutBestFitnessImprovement(n_metaepochs::Int)
        n_metaepochs > 0 || throw(ArgumentError("metaepochs number must be a positive integer, got $n_metaepochs"))
        new(n_metaepochs)
    end
end

function (lsc::MetaepochWithoutBestFitnessImprovement)(deme::Deme, metaepoch_summaries::Vector{MetaepochSummary})::Bool
    history::Vector{Float64} = deme.best_fitness_per_metaepoch
    current_best::Float64 = deme.best_fitness
    n::Int = lsc.n_metaepochs

    if length(history) < n + 1 || is_root(deme)
        return false
    end

    last_improvement_index = length(history)
    for i in length(history)-1:-1:1
        if history[i] != current_best
            break
        end
        last_improvement_index = i
    end
  
    stagnant_epochs = length(history) - last_improvement_index

    return stagnant_epochs ≥ n

end

const DEFAULT_MAX_STALE_METAEPOCH_AT_DEME = 10

DEFAULT_LSC = MetaepochWithoutBestFitnessImprovement(DEFAULT_MAX_STALE_METAEPOCH_AT_DEME)
    