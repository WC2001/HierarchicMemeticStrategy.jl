using Evolutionary
using Random

function Evolutionary.trace!(
    record::Dict{String, Any},
    objfun,
    state,
    population,
    method::Evolutionary.CMAES,
    options
)
    record["population"] = deepcopy(population)
    record["fitnesses"] = deepcopy(state.fitpop)
    record["fitness"] = deepcopy(state.fittest)
end


struct EvolutionaryCMAESMetaepoch <: MetaepochRunner end

"""
    run_metaepoch(::EvolutionaryCMAESMetaepoch, ...)

Execute the Covariance Matrix Adaptation Evolution Strategy (CMA-ES) from the `Evolutionary.jl` package.

# Supported Options
The `options::Dict{String, Any}` parameter in `TreeLevelConfig` supports the following keys when using `EvolutionaryCMAESMetaepoch`:

| Key | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `"iterations"` | `Int` | `3` | Number of CMA-ES iterations to run. |
| `"sigma0"` | `Float64` | estimate_sigma0 | Initial step size. |
| `"seed"` | `Int` | `nothing` | Local RNG seed for reproducibility. |
| `"show_trace"` | `Bool` | `false` | Whether to print `Evolutionary.jl` progress to the console. |

# Details
This method wraps the `CMAES` solver. It is particularly effective for 
"narrowing down" on an optimum once the general area
 has been found by a more exploratory algorithm like GA.

The population and fitness history are captured to populate the `MetaepochResult` trace.
"""
function run_metaepoch(
    ::EvolutionaryCMAESMetaepoch,
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

    sigma0 = get(options, "sigma0",  estimate_sigma0(initial_population))


    cmaes_params = CMAES(
        mu=Int64(population_size),  
        lambda=population_size*2,   
        sigma0=sigma0,            
    )

    cma_fitness = minimize ? (x -> fitness_function(x)) : (x -> -fitness_function(x))

    box_constraints = BoxConstraints(bounds.lower, bounds.upper)

    result = Evolutionary.optimize(
        cma_fitness,
        box_constraints,
        initial_population,
        cmaes_params,
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


function estimate_sigma0(X::Vector{Vector{Float64}})::Float64
    Xmat = reduce(vcat, [x' for x in X])
    cov_estimate = cov(Xmat; corrected=false) 
    # Compute sigma0 as sqrt(trace(cov) / dim)
    sigma0 = sqrt(tr(cov_estimate) / size(cov_estimate, 1))
    # Small epsilon to avoid zero
    eps = 1e-8
    return max(sigma0, eps)
end
