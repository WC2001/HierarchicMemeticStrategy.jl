using hms_jl: unified_population, default_sigma
using Random

@testset "Population tests" begin

    @testset "create unified population test" begin
        rng = MersenneTwister(42)
    
        lower = [-5.0, -5.0]
        upper = [5.0, 5.0]
        population_size = 10
        tree_level = 1
        sigma = [ [1.0, 1.0] for _ in 1:population_size ]

        pop = unified_population(nothing, lower, upper, population_size, tree_level, sigma, rng)

        @test length(pop) == population_size                      
        @test all(ind -> length(ind) == length(lower), pop) 
        @test all(ind -> all(i -> lower[i] ≤ ind[i] ≤ upper[i], eachindex(ind)), pop)

        lower1d = [-10.0]
        upper1d = [10.0]
        pop1d = unified_population(nothing, lower1d, upper1d, 5, tree_level, sigma[1:5], rng)
        @test all(length.(pop1d) .== 1)
        @test all(x -> lower1d[1] ≤ x[1] ≤ upper1d[1], pop1d)

        rng1 = MersenneTwister(123)
        rng2 = MersenneTwister(123)
        pop_a = unified_population(nothing, lower, upper, 5, tree_level, sigma[1:5], rng1)
        pop_b = unified_population(nothing, lower, upper, 5, tree_level, sigma[1:5], rng2)
        @test pop_a == pop_b
    end
    
end