using Evolutionary
using Random

function Evolutionary.trace!(
    record::Dict{String, Any},
    objfun,
    state,
    population,
    method::Evolutionary.GA,
    options
)
    record["population"] = deepcopy(population)
    record["fitnesses"] = deepcopy(state.fitpop)
    record["fitness"] = deepcopy(state.fitness)
end

mutable struct TraceState
    population::Vector{Vector{Vector{Float64}}}
    fitnesses::Vector{Vector{Float64}}
end

function trace_callback_factory(state_ref::Ref{TraceState})
    function callback(record)
        if record isa Vector{<:Evolutionary.OptimizationTraceRecord}
            state = state_ref[]
            last_trace = record[end]
            meta = last_trace.metadata
            if haskey(meta, "population") && haskey(meta, "fitnesses")

                push!(state.population, deepcopy(meta["population"]))
                push!(state.fitnesses, deepcopy(meta["fitnesses"]))
                    
            end
        else
            @warn "Callback received unexpected type: $(typeof(record))"
        end
        return false
    end
end

struct GaussianMutation
    sigma::Vector{Float64}
end

function (m::GaussianMutation)(individual::Vector{Float64}; rng::Random.TaskLocalRNG)
    @assert length(individual) == length(m.sigma)
    for i in eachindex(individual)
        individual[i] += randn(rng) * m.sigma[i]
    end
    return individual
end

"""
    EvolutionaryGAMetaepoch <: MetaepochRunner

A metaepoch runner that utilizes a Genetic Algorithm (GA) implementation from the 
`Evolutionary.jl` package.

# Description
This runner performs a short burst of evolution on a deme's population. It is 
typically used in the upper levels of the HMS tree for broad exploration of the 
search space.
"""
struct EvolutionaryGAMetaepoch <: MetaepochRunner end


"""
    run_metaepoch(::EvolutionaryGAMetaepoch, ...)

Execute the Genetic Algorithm for a fixed number of iterations.

# Supported Options
The `options::Dict{String, Any}` parameter in `TreeLevelConfig` supports the 
following keys when using `EvolutionaryGAMetaepoch`:

| Key | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `"iterations"` | `Int` | `3` | Number of GA iterations to run. |
| `"crossoverRate"` | `Float64` | `0.2` | Probability of crossover. |
| `"mutationRate"` | `Float64` | `0.8` | Probability of mutation. |
| `"epsilon"` | `Int` | `1` | Number of best individuals to preserve (Elitism). |
| `"crossover"` | `Function` | `LC()` | Crossover operator (from `Evolutionary.jl`). |
| `"mutation"` | `Function` | `GaussianMutation` | Custom Gaussian mutation operator. |
| `"sigma"` | `Vector` | `fill(1.0, N)` | Standard deviation for the default Gaussian mutation. |
| `"seed"` | `Int` | `nothing` | Local RNG seed for reproducibility. |
| `"show_trace"` | `Bool` | `false` | Whether to print `Evolutionary.jl` progress to console. |

# Details
The fitness function is automatically transformed based on the `minimize` flag. 
The population history is captured via a callback to ensure the `MetaepochResult` 
trace is correctly populated.
"""
function run_metaepoch(
    ::EvolutionaryGAMetaepoch,
    fitness_function::Function,
    bounds::Bounds,
    initial_population::Vector{Vector{Float64}},
    minimize::Bool,
    options::Dict{String,Any}
)::MetaepochResult

    ga_fitness = minimize ? (x->fitness_function(x)) : (x -> -fitness_function(x))
    population_size = length(initial_population)

    state_ref = Ref(TraceState(Vector{Vector{Vector{Float64}}}(), Vector{Vector{Float64}}()))
    callback = trace_callback_factory(state_ref)
    
    seed = get(options, "seed", nothing)
    rng = seed === nothing ? TaskLocalRNG() : Random.seed!(TaskLocalRNG(), seed)

    evo_opts = Evolutionary.Options(
        abstol=0.0,
        reltol=0.0,
        iterations=get(options, "iterations", 3),
        successive_f_tol=0,
        store_trace=true,
        show_trace=get(options, "show_trace", false),
        callback=callback,
        rng=rng
    )

    genome_length = length(initial_population[1])
    sigma = get(options, "sigma", fill(1.0, genome_length))

    mutation = get(options, "mutation", GaussianMutation(sigma))
    ga_params = GA(
        populationSize = population_size,
        crossoverRate = get(options, "crossoverRate", 0.2),
        mutationRate = get(options, "mutationRate", 0.8),
        epsilon = get(options, "epsilon", 1),
        crossover = get(options, "crossover", LC()),
        mutation = get(options, "mutation", mutation)
    )

    result = Evolutionary.optimize(
        ga_fitness,
        BoxConstraints(bounds.lower, bounds.upper),
        initial_population,
        ga_params,
        evo_opts
    )

    value = minimize ? result.minimum : -result.minimum
    return MetaepochResult(
        result.minimizer,
        value,
        state_ref[].population,
        state_ref[].fitnesses,
        nothing
    )
    
end

