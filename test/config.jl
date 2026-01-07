using hms_jl: default_population_sizes, default_sigma 

@testset "default population sizes" begin
    cases = [
        (1, [60]),
        (2, [60, 30]),
        (3, [60, 30, 15]),
        (4, [60, 30, 15, 8]),
        (5, [60, 30, 15, 8, 4])
    ]

    for (height, expected) in cases
        @test default_population_sizes(height) == expected
    end
end


@testset "default sigma" begin
    lower = fill(0.0, 10)
    upper = fill(100.0, 10)
    tree_height = 4

    expected = [
        fill(4.0, 10),
        fill(2.0, 10),
        fill(1.0, 10),
        fill(0.5, 10),
    ]

    result = default_sigma(lower, upper, tree_height)

    @test result == expected

end


@testset "default sigma (non-uniform domain)" begin
    lower = [0.0, -5.0, 10.0, 2.0]
    upper = [10.0,  5.0, 30.0, 6.0]
    tree_height = 3

    expected = [
        [0.4, 0.4, 0.8, 0.16],
        [0.2, 0.2, 0.4, 0.08],
        [0.1, 0.1, 0.2, 0.04],
    ]

    result = default_sigma(lower, upper, tree_height)

    @test result == expected
end