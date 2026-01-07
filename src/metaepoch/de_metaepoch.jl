using Evolutionary
using Random

function Evolutionary.trace!(
    record::Dict{String, Any},
    objfun,
    state,
    population,
    method::Evolutionary.DE,
    options
)
    record["population"] = deepcopy(population)
    record["fitnesses"] = deepcopy(state.fitness)
    record["fitness"] = deepcopy(state.fittest)
end

struct EvolutionaryDEMetaepoch <: MetaepochRunner end


"""
    run_metaepoch(::EvolutionaryDEMetaepoch, ...)

Execute the Differential Evolution (DE) algorithm implementation from the `Evolutionary.jl` package.

# Supported Options
The `options::Dict{String, Any}` parameter in `TreeLevelConfig` supports the following keys when using `EvolutionaryDEMetaepoch`:

| Key | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `"iterations"` | `Int` | `3` | Number of DE iterations to run. |
| `"F"` | `Float64` | `0.9` | Scale factor (differential weight) for mutation. |
| `"n"` | `Int` | `1` | Number of differences used in the perturbation. |
| `"selection"` | `Function` | `random` | Selection strategy (e.g., `random` or `best`). |
| `"recombination"`| `Function` | `BINX(0.5)` | Recombination operator (e.g., `BINX` or `EXP`). |
| `"K"` | `Float64` | `0.5*(F+1)` | Recombination constant. |
| `"seed"` | `Int` | `nothing` | Local RNG seed for reproducibility. |
| `"show_trace"` | `Bool` | `false` | Whether to print `Evolutionary.jl` progress to the console. |

# Details
This method wraps the `DE` solver from `Evolutionary.jl`.
"""
function run_metaepoch(
    ::EvolutionaryDEMetaepoch,
    fitness_function::Function,
    bounds::Bounds,
    initial_population::Vector{Vector{Float64}},
    minimize::Bool,
    options::Dict{String, Any}
)::MetaepochResult

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

    F = get(options, "F", 0.9)

    de_params = DE(
        populationSize=Int64(population_size),
        F=get(options, "F", 0.9),          
        n=get(options, "n", 1),
        selection=get(options, "selection", random),
        recombination=get(options, "recombination", BINX(0.5)),
        K=get(options, "K", 0.5*(F+1))
    )

    de_fitness = minimize ? (x -> fitness_function(x)) : (x -> -fitness_function(x))

    box_constraints = BoxConstraints(bounds.lower, bounds.upper)

    result = Evolutionary.optimize(
        de_fitness,
        box_constraints,
        initial_population,
        de_params,
        evo_opts
    )


    best_value = minimize ? result.minimum : -result.minimum

    return MetaepochResult(
        result.minimizer,             # Best solution
        best_value,                   # Best fitness value
        state_ref[].population,       # Population history
        state_ref[].fitnesses,        # Fitness values history
        nothing                      
    )

end
    