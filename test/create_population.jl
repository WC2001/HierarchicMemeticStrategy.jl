using Random
using HierarchicMemeticStrategy: reseed!

@testset "Population Creator Interface Tests" begin

    lower = [-5.0, -5.0]
    upper = [5.0, 5.0]
    population_size = 10
    sigma = [[1.0, 1.0], [0.5, 0.5], [0.1, 0.1]]
    seed = 42

    @testset "DefaultPopulationCreator" begin
        pc = DefaultPopulationCreator()
        reseed!(pc, MersenneTwister(seed))

        pop_root = pc(
            mean = nothing, 
            lower = lower, 
            upper = upper, 
            population_size = population_size, 
            tree_level = 1, 
            sigma = sigma
        )

        @test length(pop_root) == population_size
        @test all(ind -> all(lower .≤ ind .≤ upper), pop_root)
        
        mean_val = [0.0, 0.0]
        pop_deme = pc(
            mean = mean_val, 
            lower = lower, 
            upper = upper, 
            population_size = population_size, 
            tree_level = 2, 
            sigma = sigma
        )

        @test length(pop_deme) == population_size
        @test pop_deme[end] == mean_val
        @test all(ind -> all(lower .≤ ind .≤ upper), pop_deme)
    end

    @testset "NormalPopulationCreator" begin
        pc = NormalPopulationCreator()
        reseed!(pc, MersenneTwister(seed))
        
        mean_val = [1.0, 1.0]
        
        pop = pc(
            mean = mean_val, 
            lower = lower, 
            upper = upper, 
            population_size = population_size, 
            tree_level = 1, 
            sigma = sigma
        )

        @test length(pop) == population_size
        @test pop[end] == mean_val
    end

    @testset "Reproducibility and reseed!" begin
        pc1 = DefaultPopulationCreator()
        pc2 = DefaultPopulationCreator()
        
        reseed!(pc1, MersenneTwister(123))
        reseed!(pc2, MersenneTwister(123))

        args = (mean=nothing, lower=lower, upper=upper, population_size=5, tree_level=1, sigma=sigma)
        
        @test pc1(; args...) == pc2(; args...)
    end

    @testset "Custom Kwargs Passthrough" begin
        pc = DefaultPopulationCreator()
        @test_nowarn pc(
            mean = nothing, 
            lower = lower, 
            upper = upper, 
            population_size = 5, 
            tree_level = 1, 
            sigma = sigma,
            custom_param = "test"
        )
    end
end