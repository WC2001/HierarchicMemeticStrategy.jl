using HierarchicMemeticStrategy
using Test

@testset "HierarchicMemeticStrategy.jl" begin

    for tests in [
        "config.jl",
        "cache.jl",
        "deme.jl",
        "create_population.jl",
        "distance_metrics.jl",
        "sprout_condition.jl",
        "stop_conditions.jl",
        "evolutionary_test.jl",
        "hms.jl",
        "optimization_interface.jl"
        
    ]
        include(tests)
    end

end
