abstract type GlobalStopCondition end

"""
    GlobalMetaepochLimitReached(limit::Int)

A `GlobalStopCondition` that terminates the optimization after a fixed number of metaepochs.

# Fields
- `limit::Int`: The maximum number of metaepochs to perform (must be positive).

# Description
In HMS metaepoch represents one full iteration of the algorithm where 
every active deme in the tree is processed.

# Example
```julia
# Stop the algorithm after 50 metaepochs of the HMS tree
gsc = GlobalMetaepochLimitReached(50)
```
"""
struct GlobalMetaepochLimitReached <: GlobalStopCondition
    limit::Int

    function GlobalMetaepochLimitReached(limit::Int)
        limit > 0 || throw(ArgumentError("limit must be a positive integer, got $limit"))
        new(limit)
    end
end

function (stop_condition::GlobalMetaepochLimitReached)(metaepoch_summaries::Vector{MetaepochSummary})::Bool

    return length(metaepoch_summaries) >= stop_condition.limit
end


"""
    ProblemEvaluationLimitReached(limit::Int)

A `GlobalStopCondition` that terminates the optimization once a specific number of 
objective function evaluations has been reached.

# Fields
- `limit::Int`: The maximum allowed number of function evaluations (must be positive).

# Description
This is the most common stopping criterion for the HMS algorithm. It monitors the 
global evaluation counter and returns `true` as soon as the total evaluations 
exceed the specified `limit`, ensuring the algorithm respects a predefined 
computational budget.

# Example
```julia
# Stop the algorithm after 10,000 function evaluations
gsc = ProblemEvaluationLimitReached(10000)
```
"""
struct ProblemEvaluationLimitReached <: GlobalStopCondition
    limit::Int

    function ProblemEvaluationLimitReached(limit::Int)
        limit > 0 || throw(ArgumentError("limit must be a positive integer, got $limit"))
        new(limit)
    end
end

function (stop_condition::ProblemEvaluationLimitReached)(metaepoch_summaries::Vector{MetaepochSummary})::Bool
    return length(metaepoch_summaries) > 0 && 
            metaepoch_summaries[end].fitness_evaluation_count >= stop_condition.limit
end


const DEFAULT_MAX_FUN = 10000
GSC = ProblemEvaluationLimitReached(DEFAULT_MAX_FUN)