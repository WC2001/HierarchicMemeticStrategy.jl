using Test
using Evolutionary
using hms_jl: default_sigma, unified_population


@testset "Evolutionary GA Integration" begin

    fitness(x) = (x[1]*x[1]+x[1])*cos(x[1])

    mean = [0.0]
    lower = [-10.0]
    upper = [10.0] 
    population_size = 1000 
    tree_level = 1  
    sigma = default_sigma(lower,upper,1) 

    initial_population = unified_population(mean, lower, upper, population_size, tree_level, sigma)

    function Evolutionary.trace!(record::Dict{String, Any}, objfun, state, population, method, options)

        record["population"] = population
        
        record["fitnesses"] = state.fitpop

        record["fitness"] = state.fitness
    
    end

    options = Evolutionary.Options(
        abstol=0.0,
        reltol=0.0,
        iterations=5,
        successive_f_tol=0,
        store_trace = true,
        show_trace = false,
    )

    ga_params = GA(
        populationSize = 1000,
        crossoverRate = 0.2,
        mutationRate = 0.8,
        epsilon = 1
    )

    result = Evolutionary.optimize(
        fitness,
        BoxConstraints(lower, upper),
        initial_population,
        ga_params,
        options
    )

    
    expected_solution = [9.62]
    @test euclidean(expected_solution, result.minimizer) < 1e-1
end

