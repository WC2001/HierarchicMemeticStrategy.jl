module HierarchicMemeticStrategy
    include("cache.jl")
    include("optimization_problem.jl")
    include("individual.jl")
    include("population.jl")
    include("metaepoch/metaepoch.jl")
    include("deme/deme.jl")
    include("metaepoch/metaepoch_summary.jl")
    include("stop_conditions/stop_condition.jl")
    include("stop_conditions/global_stop_condition.jl")
    include("stop_conditions/local_stop_condition.jl")
    include("sprout_condition.jl")
    include("metaepoch/ga_metaepoch.jl")
    include("metaepoch/cmaes_metaepoch.jl")
    include("metaepoch/de_metaepoch.jl")
    include("config.jl")
    include("local_method/local_optimizer.jl")
    include("local_method/optim_lbfgs_optimizer.jl")
    include("visualization/tree_visualizer.jl")
    include("visualization/hms_result_visualizer.jl")
    include("hms_result.jl")
  
    
    using Base.Threads
    
    export HMSOptimizationProblem,
        FunctionProblem,
        hms,
        HMSResult,
        plotBestFitness,
        plotPopulations,
        plotDeme,
        summary,
        f_calls,
        iterations,
        solution,
        best_fitness,
        metaepoch_data,
        TreeLevelConfig,
        EvolutionaryGAMetaepoch,
        EvolutionaryCMAESMetaepoch,
        EvolutionaryDEMetaepoch,
        GlobalStopCondition,
        LocalStopCondition,
        ProblemEvaluationLimitReached,
        GlobalMetaepochLimitReached,
        LocalProblemEvaluationLimitReached,
        AllChildrenStopped,
        MetaepochWithoutBestFitnessImprovement,
        GaussianMutation,
        euclidean,
        manhattan,
        CustomCache,
        iscached,
        DEFAULT_LEVEL_CONFIG,
        run_metaepoch,
        Bounds,
        MetaepochRunner,
        MetaepochResult,
        LocalOptimizer,
        HMSSolver,
        AbstractPopulationCreator,
        DefaultPopulationCreator,
        NormalPopulationCreator,
        AbstractSproutCondition,
        MetricSproutCondition,
        DefaultSproutCondition

    

    """
        hms(; kwargs...)

    Execute the Hierarchical Multi-Strategy (HMS) optimization algorithm.

    # Arguments
    - `optimization_problem::HMSOptimizationProblem`: The objective function and domain bounds.
    - `level_config::Vector{TreeLevelConfig}`: Configuration for each level of the HMS tree (default: `DEFAULT_LEVEL_CONFIG`).
    - `sigma::Vector{Vector{Float64}}`: The \$ \\sigma \$ parameter for evolutionary algorithms and Gaussian mutation. It is specified for each dimension at each level of the hierarchy.
    - `gsc::GlobalStopCondition`: Condition to stop the entire optimization process.
    - `lsc::LocalStopCondition`: Condition to stop individual demes.
    - `sprout_condition::AbstractSproutCondition`: Logic determining when a deme should sprout a child deme.
    - `population_sizes::Vector{Int}`: Number of individuals in demes at each level.
    - `create_population::AbstractPopulationCreator`: Strategy object for initializing deme populations.
    - `local_optimizer::LocalOptimizer`: The specific local optimizer to use (default `LBFGSOptimizer`).
    - `use_local_method::Bool`: Whether to apply a local optimization (default: `true`).
    - `minimize::Bool`: Set to `true` to minimize, `false` to maximize (default: `true`).
    - `log_level::Int`: Verbosity of the output (0 for silent, 1 basic output, 2 basic string representation of hms tree).
    - `parallel::Bool`: Enable multi-threaded execution (default: `false`).
    - `seed::Union{Nothing, Int}`: Random seed for reproducibility.

    # Returns
    - `HMSResult`: An object containing the optimization trace, best solution, and visualization tools.

    # Example
    ```julia
    # optimizing eggholder function
    function eggholder(x::Vector{Float64})
        x1, x2 = x[1], x[2]
        return -(x2 + 47) * sin(sqrt(abs(x2 + x1/2 + 47))) - x1 * sin(sqrt(abs(x1 - (x2 + 47))))
    end

    seed = 42
    lower = [-512.0, -512.0]
    upper = [512.0, 512.0]
    problem = FunctionProblem(eggholder, lower, upper, false)
    sigma = [[100.0, 100.0], [60.0, 60.0]]
    level_config = [
        TreeLevelConfig(EvolutionaryGAMetaepoch, Dict("seed" => seed)),
        TreeLevelConfig(EvolutionaryCMAESMetaepoch, Dict("seed" => seed)),
    ]
    result = hms(
        optimization_problem = problem,
        level_config = level_config,
        sigma=sigma,
        seed=seed 
    )
    ```
    """
    function hms(;
        optimization_problem::HMSOptimizationProblem,
        level_config::Vector{TreeLevelConfig} = DEFAULT_LEVEL_CONFIG,
        sigma::Union{Nothing, Vector{Vector{Float64}}} = nothing, 
        gsc::GlobalStopCondition = GSC,
        lsc::LocalStopCondition = DEFAULT_LSC,
        sprout_condition::Union{Nothing, AbstractSproutCondition} = nothing,
        population_sizes::Union{Nothing, Vector{Int}} = nothing,
        seed::Union{Nothing, Int} = nothing,
        create_population::Union{Nothing, AbstractPopulationCreator} = nothing,
        use_local_method::Bool = true,
        local_optimizer::LocalOptimizer = LBFGSOptimizer(),
        minimize::Bool = true,
        log_level::Int = 0,
        parallel::Bool = false
    )   

        tree_height = length(level_config)

        if isnothing(sigma)
            sigma = default_sigma(optimization_problem._bounds.lower, optimization_problem._bounds.upper, tree_height)
        end

        if isnothing(population_sizes)
            population_sizes = default_population_sizes(tree_height)
        end

        if isnothing(sprout_condition)
            sprout_condition = DefaultSproutCondition(sigma)
        end

        if isnothing(create_population)
            create_population = DefaultPopulationCreator()
        end

        if !isnothing(seed)
            rng = Random.MersenneTwister(seed)
            reseed!(create_population, rng)
        end

        if tree_height < 1
            throw(ArgumentError("Tree height must be greater or equal 1."))
        end
        
        start_time = time()

        root = Deme(optimization_problem, nothing, population_sizes[1], create_population, sigma)
        demes::Vector{Deme} = [root]
        metaepochs_count::Int = 0
        metaepoch_summaries = Vector{MetaepochSummary}()
        total_evolution_time = Atomic{Float64}(0.0)
        fitness_evaluations_count = Atomic{Int}(0)

        cached_f = CustomCache(optimization_problem.fitness_function)
        
        while (!gsc(metaepoch_summaries))
            active_demes = filter(d -> d.is_active, demes)
            if length(active_demes) == 0
                println("HMS stopped as there are no active demes.")
                break
            end

            next_metaepoch_demes = filter(d -> !d.is_active, demes)
            blocked_sprouts = Vector{Vector{Float64}}()

            if parallel && Threads.nthreads() > 1
                thread_local_next = [Deme[] for _ in 1:nthreads()]
                thread_local_blocked = [Vector{Vector{Float64}}() for _ in 1:nthreads()]
                @threads for i in eachindex(active_demes)
                    tid = threadid()
                    deme = active_demes[i]

                    
                    deme_evaluations_count = 0

                    deme_f = function(x)
                        if !iscached(cached_f, x)
                            atomic_add!(fitness_evaluations_count, 1)
                            deme_evaluations_count += 1
                        end
                        return cached_f(x)
                    end

                    config = level_config[deme.level]
                    bounds = problem_bounds(optimization_problem)
                    initial_population = deepcopy(deme.population.genomes)
                    start_metaepoch_time = time()
                    metaepoch_result = run_metaepoch_at_level(
                        config,
                        deme_f,
                        bounds,
                        initial_population,
                        minimize,
                        sigma[deme.level]
                    )
                    end_metaepoch_time = time()
                    atomic_add!(total_evolution_time, end_metaepoch_time - start_metaepoch_time)

                    deme.evaluations_count += deme_evaluations_count

                    if metaepoch_result === nothing
                        deme.is_active = false
                        push!(thread_local_next[tid], deme)
                        continue
                    end

                    update_deme!(deme, metaepoch_result, minimize)

                    if lsc(deme, metaepoch_summaries)
                        deme.is_active = false
                        push!(thread_local_next[tid], deme)
                        continue
                    else
                        push!(thread_local_next[tid], deme)
                    end

                    if deme.level < tree_height
                        candidate_demes = vcat(
                            thread_local_next[tid],
                            get_not_processed_demes(demes, thread_local_next[tid])
                        )

                        if sprout_condition(metaepoch_result.solution, deme.level + 1, candidate_demes)
                            new_deme = Deme(
                                optimization_problem,
                                deme,
                                population_sizes[deme.level + 1],
                                create_population,
                                sigma
                            )
                            push!(thread_local_next[tid], new_deme)
                        else
                            push!(thread_local_blocked[tid], metaepoch_result.solution)
                        end
                    end
                end
                next_metaepoch_demes = vcat(next_metaepoch_demes, reduce(vcat, thread_local_next))
                blocked_sprouts = reduce(vcat, thread_local_blocked)
        
            else
            
                for deme in active_demes

                    deme_evaluations_count = 0
        
                    deme_f = function(x)
                        
                        if !iscached(cached_f, x)
                            atomic_add!(fitness_evaluations_count, 1)
                            deme_evaluations_count += 1
                        end
                        return cached_f(x)
                    end
            
                    config::TreeLevelConfig = level_config[deme.level]
                    bounds::Bounds = problem_bounds(optimization_problem)
                    initial_population::Vector{Vector{Float64}} = deepcopy(deme.population.genomes)
                    start_metaepoch_time = time()
                    metaepoch_result::MetaepochResult = run_metaepoch_at_level(
                        config,
                        deme_f,
                        bounds,
                        initial_population, 
                        minimize,
                        sigma[deme.level]
                    )
                    end_metaepoch_time = time()
                    atomic_add!(total_evolution_time, end_metaepoch_time - start_metaepoch_time)
        
                    deme.evaluations_count += deme_evaluations_count
        
                    if metaepoch_result === nothing
                        deme.is_active = false
                        push!(next_metaepoch_demes, deme)
                        continue
                    end
        
                    deme = update_deme!(deme, metaepoch_result, minimize)
                    
                    if lsc(deme, metaepoch_summaries)
                        deme.is_active = false
                        push!(next_metaepoch_demes, deme)
                        continue
                    else
                        push!(next_metaepoch_demes, deme)
                    end
        
                    if deme.level >= tree_height
                        continue
                    end
        
                    candidate_demes = vcat(next_metaepoch_demes, get_not_processed_demes(demes, next_metaepoch_demes))
        
                    if sprout_condition(metaepoch_result.solution, deme.level + 1, candidate_demes)
                        new_deme = Deme(optimization_problem, deme, population_sizes[deme.level + 1], create_population, sigma)
                        push!(next_metaepoch_demes, new_deme)
                    else
                        push!(blocked_sprouts, metaepoch_result.solution)
                    end
                end
            end

            demes = next_metaepoch_demes
            best = get_best_solution(demes, minimize)
            summary = MetaepochSummary(
                deepcopy(demes),
                best.fitness,
                best.solution,
                total_evolution_time[],
                fitness_evaluations_count[],
                blocked_sprouts,
            )

            if log_level > 0
                log_metaepoch_summary(summary, metaepochs_count, root, best, log_level > 1)
            end
            push!(metaepoch_summaries, summary)
            metaepochs_count += 1

        end

        if use_local_method
            log_level > 0  && println("Starting local method.")

            leaves::Vector{Deme} = get_leaves(demes, tree_height)
            
            for leaf in leaves
                local_method_f_calls = optimize!(local_optimizer, leaf, cached_f)
                atomic_add!(fitness_evaluations_count, local_method_f_calls)
            end

            best = get_best_solution(demes, minimize)

            summary = MetaepochSummary(
                demes,
                best.fitness,
                best.solution,
                0.0,
                fitness_evaluations_count[],
                [],
            )

            if log_level > 0
                log_metaepoch_summary(summary, metaepochs_count, root, best, log_level > 1)
            end
            push!(metaepoch_summaries, summary)
            metaepochs_count += 1
        end

        result_visualizer = HMSResultVisualizer(metaepoch_summaries)
        total_time = time() - start_time

        hms_result = HMSResult(
            metaepoch_summaries,
            result_visualizer,
            optimization_problem
        )
        
        return hms_result
    end

    using Optimization
    using OptimizationBase: allowsbounds
    using SciMLBase: ReturnCode.Success

    """
        HMSSolver(; kwargs...)

    Optimization algorithm for the `Optimization.jl` interface.

    # Parameters
    - `level_config::Vector{TreeLevelConfig}`: Configuration for each level of the HMS tree.
    - `sigma::Vector{Vector{Float64}}`: The \$ \\sigma \$ parameter for evolutionary algorithms and Gaussian mutation. It is specified for each dimension at each level of the hierarchy.
    - `gsc::GlobalStopCondition`: Condition to stop the entire optimization process.
    - `lsc::LocalStopCondition`: Condition to stop individual demes.
    - `sprout_condition::AbstractSproutCondition`: Logic determining when a deme should sprout a child deme.
    - `population_sizes::Vector{Int}`: Number of individuals in demes at each level.
    - `create_population::AbstractPopulationCreator`: Strategy object for initializing deme populations.
    - `local_optimizer::LocalOptimizer`: The specific local optimizer to use (default `LBFGSOptimizer`).
    - `use_local_method::Bool`: Whether to apply a local optimization (default: `true`).
    - `minimize::Bool`: Set to `true` to minimize, `false` to maximize (default: `true`).
    - `log_level::Int`: Verbosity of the output (0 for silent, 1 basic output, 2 basic string representation of hms tree).
    - `parallel::Bool`: Enable multi-threaded execution (default: `false`).
    - `seed::Union{Nothing, Int}`: Random seed for reproducibility.
    """
    Base.@kwdef struct HMSSolver <: Optimization.SciMLBase.AbstractOptimizationAlgorithm
        level_config::Vector{TreeLevelConfig} = DEFAULT_LEVEL_CONFIG
        sigma::Union{Nothing, Vector{Vector{Float64}}} = nothing
        gsc::GlobalStopCondition = GSC
        lsc::LocalStopCondition = DEFAULT_LSC
        sprout_condition::Union{Nothing, AbstractSproutCondition} = nothing
        population_sizes::Union{Nothing, Vector{Int}} = nothing
        seed::Union{Nothing, Int} = nothing
        create_population::Union{Nothing, AbstractPopulationCreator} = nothing
        use_local_method::Bool = true
        local_optimizer::LocalOptimizer = LBFGSOptimizer()
        log_level::Int = 0
        parallel::Bool = false
        minimize::Bool = true
    end

    OptimizationBase.allowsbounds(::HMSSolver) = true

    function Optimization.SciMLBase.__solve(prob::Optimization.SciMLBase.OptimizationProblem, opt::HMSSolver, args...; 
                                       kwargs...)

        if isnothing(prob.lb) || isnothing(prob.ub)
            error("HMSSolver requires a bounded optimization problem (lb and ub must be provided).")
        end
        
        cache = SciMLBase.init(prob, opt, args...; kwargs...)

        hms_prob = FunctionProblem(
            fitness_function = (x; kwargs...) -> prob.f(x, prob.p),
            lower = prob.lb,
            upper = prob.ub,
            maximize = !opt.minimize,
            u0 = prob.u0
        )

        hms_res = hms(
            optimization_problem = hms_prob,
            level_config = opt.level_config,
            sigma = opt.sigma,
            gsc = opt.gsc,
            lsc = opt.lsc,
            sprout_condition = opt.sprout_condition,
            population_sizes = opt.population_sizes,
            seed = opt.seed,
            create_population = opt.create_population,
            use_local_method = opt.use_local_method,
            local_optimizer = opt.local_optimizer,
            minimize = opt.minimize,
            log_level = opt.log_level,
            parallel = opt.parallel
        )

        return SciMLBase.build_solution(
            cache, 
            opt, 
            hms_res.solution, 
            hms_res.best_fitness;
            original = hms_res,
            retcode = SciMLBase.ReturnCode.Success
        )

    end

end
