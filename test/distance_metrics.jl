@testset "Distance metrics" begin
    a = [1.0, 2.0, 3.0]
    b = [4.0, 6.0, 3.0]

    # Euclidean
    @test euclidean(a, a) == 0.0
    @test euclidean(b, b) == 0.0
    @test euclidean(a, b) == sqrt(25)

    # Manhattan
    @test manhattan(a, a) == 0.0
    @test manhattan(b, b) == 0.0
    @test manhattan(a, b) == 7.0 
end