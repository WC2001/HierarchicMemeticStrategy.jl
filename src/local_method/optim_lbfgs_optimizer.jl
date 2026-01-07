using Optim
using LineSearches

Base.@kwdef struct LBFGSOptimizer <: LocalOptimizer
    max_iters::Int = 5
    options::Dict{String,Any} = Dict{String,Any}()
end

function optimize!(opt::LBFGSOptimizer, deme::Deme, cached_f::CustomCache)

    problem::OptimizationProblem = deme.population.problem
    deme_evaluations_count = 0
    deme_f = function(x)
        if !iscached(cached_f, x)
            deme_evaluations_count += 1
        end
        return cached_f(x)
    end
    f = maximize_fitness(problem) ? (x->-deme_f(x)) : (x ->deme_f(x))
    x0 = deme.best_solution
    bounds = problem_bounds(problem)

    lbfgs_kwargs = Dict{Symbol, Any}()
    for (k, v) in opt.options
        lbfgs_kwargs[Symbol(k)] = v
    end
    

    if !haskey(lbfgs_kwargs, :linesearch)
        hz = LineSearches.HagerZhang(
            delta=0.1,
            sigma=0.9,
            alphamax=Inf,
            rho=5.0,
            epsilon=1e-6,
            gamma=0.66,
            linesearchmax=50,
            psi3=0.1,
            display=0
        )

        lbfgs_kwargs[:linesearch] = hz
    end

    kwargs_namedtuple = (; lbfgs_kwargs...)

    lbfgs = LBFGS(; kwargs_namedtuple...)
    fminbox = Fminbox(lbfgs)

    result = optimize(
        f,
        bounds.lower,
        bounds.upper,
        x0,
        fminbox,
        Optim.Options(
            outer_iterations = 5,
            iterations = 10,
            x_abstol = 1e-8,     
            f_reltol = 1e-8,
            f_calls_limit = 10,
            outer_x_abstol = 1e-8,
            outer_f_reltol = 1e-8
        )
    )

    # println("iterations: ", Optim.iterations(result))
    # println(Optim.iteration_limit_reached(result))
    # println("converged: ", Optim.converged(result))
    # println("f_calls: ", Optim.f_calls(result), " g_calls: ", Optim.g_calls(result))
    # println("fitness calls: ", deme_evaluations_count)
   
    deme.evaluations_count += deme_evaluations_count
    update!(deme, result)
       
    return deme_evaluations_count
end

