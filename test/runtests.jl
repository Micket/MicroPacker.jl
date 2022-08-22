using MicroPacker
using Test
using StableRNGs

#u, u₀ = Uniform(0.25, 0.75), 0.5
#w, w₀ = Weibull(0.1, 0.5), 0.5
c, c₀ = Dirac(0.5), 0.5

const ttgenerator = TTGenerator(c, c, c)
const cgenerator = CylinderGenerator(c, c)
const sgenerator = SphereGenerator(c)

@testset "GenerateShape" begin
    @test 1 == 1 # todo
end

@testset "Volume" begin
    @test 1 == 1 # todo
end

@testset "Inside" begin
    @test 1 == 1 # todo
end