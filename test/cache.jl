using Base.Threads: Atomic, atomic_add!
using hms_jl: CustomCache, iscached

@testset "CustomCache" begin

    rosenbrock(x) = (1 - x[1])^2 + (100 * (x[2] - x[1]^2)^2)

    cached_rosenbrock = CustomCache(rosenbrock)
    input1 = [1.0, 1.0]
    input2 = [1.0, 1.1]
    input3 = [1.0, 1.0]
    input4 = [1.0, 1.1]

    val1 = cached_rosenbrock(input1)
    val2 = cached_rosenbrock(input2)
    @test iscached(cached_rosenbrock, input3)
    @test val1 == cached_rosenbrock(input3)
    @test iscached(cached_rosenbrock, input4)
    @test val2 == cached_rosenbrock(input4)


end

@testset "CustomCache Concurrent with Call Counter" begin
    callcount = Atomic{Int}(0)

    function rosenbrock(x)
        atomic_add!(callcount, 1)
        return (1 - x[1])^2 + (100 * (x[2] - x[1]^2)^2)
    end

    cached_rosenbrock = CustomCache(rosenbrock)

    inputs = [[1.0, 1.0], [1.0, 1.1], [0.9, 1.2]]

    expected = Dict(x => rosenbrock(x) for x in inputs)

    callcount[] = 0

    nthreads = Threads.nthreads()
    niter = 200

    @sync for tid in 1:nthreads
        Threads.@spawn begin
            for i in 1:niter
                x = inputs[rand(1:length(inputs))]
                v = cached_rosenbrock(x)

                @test v == expected[x]
                @test iscached(cached_rosenbrock, x)
            end
        end
    end

    @test callcount[] <= length(inputs)
end