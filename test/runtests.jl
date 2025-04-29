using TrixiShockTube
using Test
using Unitful
using PyThermo
using PyThermo.ShockTube

@testset "TrixiShockTube.jl" begin
    @testset "docs/src/index.md example replication" begin
        N₂ = Species("N2"; P = 3.0u"bar")
        N₂Ar = Mixture(["N2" => 0.5, "Ar" => 0.5])
        @test pressure(N₂Ar) == 101325u"Pa"

        slabs = (N₂   => ustrip(u"m", 30.0u"inch"),
            N₂Ar => ustrip(u"m", 125.0u"inch"))
        ic, equations = TrixiShockTube.build_shocktube_ic(slabs)
        semi = build_semidiscretization(slabs, ic, equations;
            initial_refinement_level=8)
        saveat = LinRange(0, 0.02, 1000)
        limiter! = build_limiter(equations)
        sol = run_shock(semi, saveat, 0.5, limiter!);
        x, t, data = xtdata(sol, semi)

        Mₛ = calculate_Mach(N₂, N₂Ar)
        sc = shockcalc(N₂, N₂Ar, Mₛ)
        @test isapprox(density(sc.reflected), 
            maximum(data[:rho2][end, :])*u"kg/m^3", rtol=5e-2)
    end
end
