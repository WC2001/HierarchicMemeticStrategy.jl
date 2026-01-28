using Optimization

@testset "Optimization.jl Interface Integration" begin
    seed = 42

    @testset "Basic Quadratic (1D)" begin
        f(x, p) = x[1]*x[1] + x[1] + 2
        lb = [-5.0]
        ub = [5.0]
        prob = OptimizationProblem(f, [0.0], lb=lb, ub=ub)
        
        sol = solve(prob, HMSSolver(seed=seed))
        
        expected_result = [-0.5]
        @test euclidean(sol.u, expected_result) < 1e-2
        @test sol.original isa HMSResult
    end

    @testset "Rosenbrock 2D" begin
        f(x, p) = (p[1] - x[1])^2 + (p[2] * (x[2] - x[1]^2)^2)
        lb = [-30.0, -30.0]
        ub = [30.0, 30.0]
        p = [1.0, 100.0]
        
        level_config = [
            TreeLevelConfig(EvolutionaryGAMetaepoch, Dict("seed" => seed)),
            TreeLevelConfig(EvolutionaryCMAESMetaepoch, Dict("seed" => seed)),
        ]
        
        prob = OptimizationProblem(f, [0.0, 0.0], p, lb=lb, ub=ub)
        sol = solve(prob, HMSSolver(level_config=level_config, seed=seed))
        @test euclidean(sol.u, [1.0, 1.0]) < 1e-2
    end

    @testset "Ackley 10D" begin
        function ackley(x, p)
            a, b, c = 20, 0.2, 2π
            d = length(x)
            sum_sq = sum(xi^2 for xi in x)
            sum_cos = sum(cos(c * xi) for xi in x)
            return -a * exp(-b * sqrt(sum_sq / d)) - exp(sum_cos / d) + a + exp(1)
        end
    
        D = 10
        lb = fill(-32.768, D)
        ub = fill(32.768, D)
        prob = OptimizationProblem(ackley, fill(0.0, D), lb=lb, ub=ub)
    
        sol = solve(prob, HMSSolver(seed=seed))
        
        @test length(sol.u) == 10
        @test euclidean(fill(0.0, D), sol.u) < 1.0
    end

    @testset "Eggholder" begin
        function eggholder(x, p)
            x1, x2 = x[1], x[2]
            return -(x2 + 47) * sin(sqrt(abs(x2 + x1/2 + 47))) - x1 * sin(sqrt(abs(x1 - (x2 + 47))))
        end

        lb = [-512.0, -512.0]
        ub = [512.0, 512.0]
        prob = OptimizationProblem(eggholder, [0.0, 0.0], lb=lb, ub=ub)
        
        sol = solve(prob, HMSSolver(seed=seed))
        
        expected_solution = [512.0, 404.2319]
        @test abs(sol.objective - eggholder(expected_solution, nothing)) < 1e2
    end

    @testset "Beale" begin
        function beale(x, p)
            return (1.5 - x[1] + x[1]*x[2])^2 +
                   (2.25 - x[1] + x[1]*x[2]^2)^2 +
                   (2.625 - x[1] + x[1]*x[2]^3)^2
        end

        lb, ub = [-4.5, -4.5], [4.5, 4.5]
        prob = OptimizationProblem(beale, [0.0, 0.0], lb=lb, ub=ub)
        sol = solve(prob, HMSSolver(seed=seed))

        @test euclidean(sol.u, [3.0, 0.5]) < 1e-4 
    end
end