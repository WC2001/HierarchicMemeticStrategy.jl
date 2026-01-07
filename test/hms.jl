@testset "hms test" begin

    seed = 42

    @testset "basic example" begin
        f = FunctionProblem(x -> x[1]*x[1]+x[1]+2, [-5], [5], false)
        result = hms(optimization_problem=f, seed=seed)

        expected_result = [-0.5]

        @test euclidean(result.solution, expected_result) < 1e-2

    end

    @testset "basic example parallel" begin
        f = FunctionProblem(x -> x[1]*x[1]+x[1]+2, [-5], [5], false)
        result = hms(optimization_problem=f, seed=seed, parallel=true)

        expected_result = [-0.5]

        @test euclidean(result.solution, expected_result) < 1e-2

    end
 
    
    @testset "rosenbrock 2D" begin
        rosenbrock(x) = (1 - x[1])^2 + (100 * (x[2] - x[1]^2)^2)
        lower = [Float64(-30), Float64(-30)]
        upper = [Float64(30), Float64(30)]
        
        level_config = [
            TreeLevelConfig(EvolutionaryGAMetaepoch, Dict("seed" => seed)),
            TreeLevelConfig(EvolutionaryCMAESMetaepoch, Dict("seed" => seed)),
        ]
        
        problem = FunctionProblem(rosenbrock, lower, upper, false)
        result = hms(optimization_problem=problem, level_config=level_config, seed=seed)

        expected_result = [1.0, 1.0]
        @test euclidean(result.solution, expected_result) < 1e-2

    end

    @testset "rosenbrock 2D v2" begin
        rosenbrock(x) = (1 - x[1])^2 + (100 * (x[2] - x[1]^2)^2)
        lower = [Float64(0), Float64(0)]
        upper = [Float64(2), Float64(2)]

        level_config = [
            TreeLevelConfig(EvolutionaryGAMetaepoch, Dict("seed" => seed)),
            TreeLevelConfig(EvolutionaryDEMetaepoch, Dict("iterations" => 3, "seed" => seed)),
        ]
        
        problem = FunctionProblem(rosenbrock, lower, upper, false)
        result = hms(optimization_problem=problem, level_config=level_config, seed=seed)
        expected_result = [1.0, 1.0]
        @test euclidean(result.solution, expected_result) < 1e-2
    end

    @testset "eggholder" begin
        function eggholder(x::Vector{Float64})
            x1, x2 = x[1], x[2]
            return -(x2 + 47) * sin(sqrt(abs(x2 + x1/2 + 47))) - x1 * sin(sqrt(abs(x1 - (x2 + 47))))
        end

        lower = [Float64(-512), Float64(-512)]
        upper = [Float64(512), Float64(512)]
        #sigma = [[200.0, 200.0], [100.0, 100.0]]
        sigma = [[100.0, 100.0], [60.0, 60.0]]
        level_config = [
            TreeLevelConfig(EvolutionaryGAMetaepoch, Dict("seed" => seed)),
            TreeLevelConfig(EvolutionaryCMAESMetaepoch, Dict("seed" => seed)),
        ]
        problem = FunctionProblem(eggholder, lower, upper, false)
        result = hms(optimization_problem=problem, level_config=level_config, sigma=sigma, seed=seed)
        expected_solution = [512, 404.2319]
        expected_fitness = eggholder(expected_solution)

        @test abs(result.best_fitness - expected_fitness) < 1e2

    end

    @testset "rastrigin" begin
        function rastrigin(x::Vector{Float64})
            A = 10
            n = length(x)
            return A * n + sum(xi^2 - A * cos(2 * π * xi) for xi in x)
        end
    
        lower = [Float64(-5.12), Float64(-5.12)]
        upper = [Float64(5.12), Float64(5.12)]
    
        level_config = [
            TreeLevelConfig(EvolutionaryGAMetaepoch, Dict("seed" => seed)),
            TreeLevelConfig(EvolutionaryCMAESMetaepoch, Dict("seed" => seed)),
        ]
    
        problem = FunctionProblem(rastrigin, lower, upper, false)
        result = hms(optimization_problem=problem, level_config=level_config, seed=seed)

        expected_solution = [0.0, 0.0]
        
        @test length(result.solution) == 2
        @test euclidean(expected_solution, result.solution) < 1e-3
        
    end

    @testset "ackley 2D" begin
        function ackley(x::Vector{Float64})
            a = 20
            b = 0.2
            c = 2π
            d = length(x)
            sum_sq = sum(xi^2 for xi in x)
            sum_cos = sum(cos(c * xi) for xi in x)
            term1 = -a * exp(-b * sqrt(sum_sq / d))
            term2 = -exp(sum_cos / d)
            return term1 + term2 + a + exp(1)
        end
    
        lower = [Float64(-32.768), Float64(-32.768)]
        upper = [Float64(32.768), Float64(32.768)]
    
        level_config = [
            TreeLevelConfig(EvolutionaryGAMetaepoch, Dict("seed" => seed)),
            TreeLevelConfig(EvolutionaryCMAESMetaepoch, Dict("seed" => seed)),
        ]
    
        problem = FunctionProblem(ackley, lower, upper, false)
        result = hms(optimization_problem=problem, level_config=level_config, seed=seed)
        
        expected_solution = [0.0, 0.0]
        
        @test length(result.solution) == 2
        @test euclidean(expected_solution, result.solution) < 1e-2
    end

    @testset "ackley 10D" begin
        function ackley(x::Vector{Float64})
            a = 20
            b = 0.2
            c = 2π
            d = length(x)
            sum_sq = sum(xi^2 for xi in x)
            sum_cos = sum(cos(c * xi) for xi in x)
            term1 = -a * exp(-b * sqrt(sum_sq / d))
            term2 = -exp(sum_cos / d)
            return term1 + term2 + a + exp(1)
        end
    
        D = 10
        lower = fill(-32.768, D)
        upper = fill(32.768, D)
    
        level_config = [
            TreeLevelConfig(EvolutionaryGAMetaepoch, Dict("seed" => seed)),
            TreeLevelConfig(EvolutionaryCMAESMetaepoch, Dict("seed" => seed)),
        ]
    
        problem = FunctionProblem(ackley, lower, upper, false)
        result = hms(optimization_problem=problem, level_config=level_config, seed=seed)
        
        expected_solution = fill(0.0, D)
        
        @test length(result.solution) == 10
        @test euclidean(expected_solution, result.solution) < 1e0

    end

    @testset "Schwefel 2D" begin
        function schwefel(x::Vector{Float64})
            n = length(x)
            418.9829 * n - sum(xi * sin(sqrt(abs(xi))) for xi in x)
        end

        lower = [-500.0, -500.0]
        upper = [500.0, 500.0]

        level_config = [
            TreeLevelConfig(EvolutionaryGAMetaepoch, Dict("seed" => seed)),
            TreeLevelConfig(EvolutionaryCMAESMetaepoch, Dict("seed" => seed)),
        ]

        problem = FunctionProblem(schwefel, lower, upper, false)
        result = hms(optimization_problem=problem, level_config=level_config, seed=seed)

        expected_solution = [420.9687, 420.9687]

        @test length(result.solution) == 2
        @test euclidean(expected_solution, result.solution) < 1.0
    end


    @testset "Griewank 2D" begin
        function griewank(x::Vector{Float64})
            n = length(x)
            sum_sq = sum(xi^2 for xi in x) / 4000
            prod_cos = prod(cos(x[i] / sqrt(i)) for i in 1:n)
            return sum_sq - prod_cos + 1
        end
        seed = 10
        lower = [-600.0, -600.0]
        upper = [600.0, 600.0]

        level_config = [
            TreeLevelConfig(EvolutionaryGAMetaepoch, Dict("seed" => seed)),
            TreeLevelConfig(EvolutionaryCMAESMetaepoch, Dict("seed" => seed)),
        ]

        problem = FunctionProblem(griewank, lower, upper, false)
        result = hms(optimization_problem=problem, level_config=level_config, seed=seed)

        expected_solution = [0.0, 0.0]

        @test length(result.solution) == 2
        @test euclidean(expected_solution, result.solution) < 1e1
    end

    @testset "Beale function" begin
        function beale(x::Vector{Float64})
            return (1.5 - x[1] + x[1]*x[2])^2 +
                (2.25 - x[1] + x[1]*x[2]^2)^2 +
                (2.625 - x[1] + x[1]*x[2]^3)^2
        end

        lower = [-4.5, -4.5]
        upper = [4.5, 4.5]

        level_config = [
            TreeLevelConfig(EvolutionaryGAMetaepoch, Dict("seed" => seed)),
            TreeLevelConfig(EvolutionaryCMAESMetaepoch, Dict("seed" => seed)),
        ]

        problem = FunctionProblem(beale, lower, upper, false)
        result = hms(optimization_problem = problem, level_config = level_config, seed = seed)

        expected_solution = [3.0, 0.5]

        @test length(result.solution) == 2
        @test euclidean(expected_solution, result.solution) < 1e-4 
    end


end