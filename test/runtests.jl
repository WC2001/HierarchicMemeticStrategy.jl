using hms_jl
using Test

@testset "hms_jl.jl" begin

    for tests in [
        "config.jl",
        "cache.jl",
        "deme.jl",
        "create_population.jl",
        "distance_metrics.jl",
        "sprout_condition.jl",
        "stop_conditions.jl",
        "evolutionary_test.jl",
        "hms.jl"
        
    ]
        include(tests)
    end

end
