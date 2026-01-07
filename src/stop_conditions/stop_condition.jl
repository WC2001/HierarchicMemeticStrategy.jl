abstract type StopCondition end

struct MetaepochLimit <: StopCondition
    limit::Int
end

function (sc::MetaepochLimit)(obj::Union{Deme})::Bool
    return length(obj.best_solution_per_metaepoch) >= sc.limit
end

struct DontStop <: StopCondition end
(sc::DontStop)(::Union{Deme})::Bool = false

struct DontRun <: StopCondition end
(sc::DontRun)(::Union{Deme})::Bool = true