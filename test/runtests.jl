using Test
using SeisSolXDMF
using Pkg.Artifacts

@testset "XDMF" begin
    test_resources_path = artifact"test_resources"

    filename = joinpath(test_resources_path, "output-surface-h5.xdmf")
    xdmf = XDMFFile(filename)

    @testset "Header" begin
        timesteps = timesteps_of(xdmf)
        @test length(timesteps) == 3
        @test timesteps == [float(i) / 10 for i ∈ 0:2]
        @test eltype(timesteps) == Float64
    end

    @testset "HDF5 - Surface" begin
        simplices, points = grid_of(xdmf)
        @test size(simplices) == (3, 172616)
        @test size(points) == (3, 86468)

        u = data_of(xdmf, 1, "u")
        @test ndims(u) == 1
        @test length(u) == size(simplices)[2]

        u = data_of(xdmf, 3, "u")
        @test minimum(u) < 0.
        @test maximum(u) > 0.
    end

    filename = joinpath(test_resources_path, "output-volume-bin.xdmf")
    xdmf = XDMFFile(filename)

    @testset "POSIX - Volume" begin
        simplices, points = grid_of(xdmf)
        @test size(simplices) == (4, 10000)
        @test size(points) == (3, 8000)

        u = data_of(xdmf, 1, "u")
        @test ndims(u) == 1
        @test length(u) == size(simplices)[2]
        @test all(tup -> begin (i, val) = tup; val == i * .001 end, enumerate(u))

        u = data_of(xdmf, 3, "u")
        @test all(tup -> begin (i, val) = tup; val ≈ 200 + i * .001 end, enumerate(u))
    end
end