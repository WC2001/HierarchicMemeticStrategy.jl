using HierarchicMemeticStrategy: default_sigma, default_create_population, Deme, MetaepochResult

@testset "Deme struct" begin

    @testset "Create 1D Deme" begin
        lower = [-10.0]
        upper = [10.0]
        tree_height = 2
        sigma = default_sigma(lower, upper, tree_height)
        population_size = 20
        create_population = default_create_population(sigma)
        f(x) = (x) -> x^3 + 2x^2 - 11 
        problem = FunctionProblem(f, lower, upper, false)
        deme = Deme(problem, nothing, population_size, create_population)

        @test deme.level == 1
        @test length(deme.population.genomes) == population_size
        @test length(deme.population.genomes[1]) == 1
        @test deme.is_active == true
        @test deme.parent_id === nothing
        @test deme.sprout === nothing
    end

    @testset "Create 2D Deme" begin
        lower = [-5.0, -5.0]
        upper = [5.0, 5.0]
        tree_height = 3
        population_size = 35
        sigma = default_sigma(lower, upper, tree_height)
        create_population = default_create_population(sigma)
        rosenbrock(x) = (1 - x[1])^2 + (100 * (x[2] - x[1]^2)^2)
        problem = FunctionProblem(rosenbrock, lower, upper, false)
       
        deme = Deme(problem, nothing, population_size, create_population)

        @test deme.level == 1
        @test length(deme.population.genomes) == population_size
        @test length(deme.population.genomes[1]) == 2
        @test deme.is_active
        @test deme.parent_id === nothing
        @test deme.sprout === nothing
    end

    @testset "Create 2D Deme with parent" begin
        lower = [-5.0, -5.0]
        upper = [5.0, 5.0]
        tree_height = 2
        population_size = 45
        sigma = default_sigma(lower, upper, tree_height)
        create_population = default_create_population(sigma)
        rosenbrock(x) = (1 - x[1])^2 + (100 * (x[2] - x[1]^2)^2)
        problem = FunctionProblem(rosenbrock, lower, upper, false)
        parent = Deme(problem, nothing, population_size, create_population)
        sprout = [1.5, -2.2]
        parent.best_solution_per_metaepoch = [sprout]
        child = Deme(problem, parent, population_size, create_population)

        @test child.level == 2
        @test length(child.population.genomes) == population_size
        @test length(child.population.genomes[1]) == 2
        @test child.parent_id == parent.id
        @test child.sprout == sprout
    end


    @testset "Update Deme" begin
        lower = [-5.0, -5.0]
        upper = [5.0, 5.0]
        population_size = 10

        sigma = default_sigma(lower, upper, 2)
        deme = Deme()
       
        solution = [1.0, 2.0]
        value = 42.0
        population = [solution for _ in 1:population_size]
        fitnesses = [value for _ in 1:population_size]
        metaepoch_result = MetaepochResult(
            solution,
            42.0,
            [population],
            [fitnesses],
            nothing
        )

        HierarchicMemeticStrategy.update_deme!(deme, metaepoch_result, false)

        @test deme.best_solution == solution
        @test deme.best_fitness == value
        @test deme.population.genomes == population
        @test deme.best_fitness_per_metaepoch == [value]
        @test deme.best_solution_per_metaepoch == [solution]
        @test deme.is_active
    end
end