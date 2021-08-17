using Test
using SeisSolXDMF

@testset "XDMF" begin
    workdir = pwd()
    if endswith(workdir, "test")
        cd("..")
    end

    filename = "test_resources/output-surface.xdmf"
    xdmf = XDMFFile(filename)
    timesteps = timesteps_of(xdmf)
    @test length(timesteps) == 3
    @test timesteps == [float(i) / 10 for i ∈ 0:2]
    @test eltype(timesteps) == Float64

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