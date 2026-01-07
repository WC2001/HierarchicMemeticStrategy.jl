using hms_jl: default_sigma, sprout_default_euclidean_distances

@testset "default sprout condition distances" begin
    lower = fill(0.0, 10)
    upper = fill(100.0, 10)
    tree_height = 3

    sigma = default_sigma(lower, upper, tree_height)
    distances = sprout_default_euclidean_distances(sigma)

    expected = [24.0, 12.0, 6.0]

    @test distances ≈ expected
end


@testset "default sprout condition distances (non-uniform domain)" begin
    lower = [0.0, -5.0, 10.0, 2.0]
    upper = [10.0,  5.0, 30.0, 6.0]
    tree_height = 3

    sigma = default_sigma(lower, upper, tree_height)
    distances = sprout_default_euclidean_distances(sigma)

    expected = [1.056, 0.528, 0.264]

    @test distances ≈ expected
end