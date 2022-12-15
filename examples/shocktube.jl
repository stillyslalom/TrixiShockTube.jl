### A Pluto.jl notebook ###
# v0.19.18

using Markdown
using InteractiveUtils

# ╔═╡ f1b8cf12-7bd3-11ed-0e4d-6934b16fe50f
begin
	using Pkg
	Pkg.activate(".")
	using Revise
	using TrixiShockTube
	using TrixiShockTube.Trixi: varnames, PlotData1D, cons2prim, LobattoLegendreBasis
	using PlutoUI
	using CairoMakie
	using Unitful
	using PyThermo
	using PyThermo.ShockTube
end

# ╔═╡ b9cc394e-5951-40ae-b53d-f9454fe20421
driver = Mixture(["SF6" => 0.93, "acetone" => 0.07], P = 45u"psi", T=22u"°C")

# ╔═╡ f3c67e3f-6b14-4493-a9e1-71754375c30b
driven = Mixture(["SF6" => 0.93, "acetone" => 0.07], P = 14.5u"psi", T=22u"°C")

# ╔═╡ 3cee7d89-132e-4b15-845e-e9ec510539e3
ambient = Species("N2", P = 14.5u"psi", T=19u"°C")

# ╔═╡ 68b78569-9c4a-4b79-b6cd-089aa53e0de1
Mₛ = calculate_Mach(driver, driven)

# ╔═╡ edc0fb1c-b9e4-4274-bf08-aba90188372b
sc = shockcalc(driver, driven, Mₛ)

# ╔═╡ 756d3adb-b79d-4737-9e18-0edf08b9ad5f
slabs = (driver => ustrip(u"m", 3.308u"inch"),
		 driven => ustrip(u"m", 10.61u"inch"), 
		 ambient => ustrip(u"m", 30u"inch"))

# ╔═╡ b8d1b888-d501-47a8-a0b3-670ea10ad6b0
sum(last.(slabs)[1:2])

# ╔═╡ d87dec7b-7f0e-402e-900b-a9263c894a87
TrixiShockTube.build_shocktube_ic

# ╔═╡ 10bda30d-a84f-41d9-8c21-585f976cbcb7
ic, equations = TrixiShockTube.build_shocktube_ic(slabs)

# ╔═╡ 14b66293-209b-4e61-b7dc-f52755bf3378
semi = build_semidiscretization(slabs, ic, equations;
	basis = LobattoLegendreBasis(3), initial_refinement_level=8)

# ╔═╡ d0126bd1-24c3-45f4-97e7-39a4fc4e00ae
limiter! = build_limiter(equations)

# ╔═╡ 6b0c4c1d-f81e-4357-ba24-679c2c06ded9
ode = build_ODE(semi, LinRange(0, 0.04, 1024), 0.5; callbacks=())

# ╔═╡ 95ff2fa9-d64b-4aa2-bee2-e708c64edd70
sol = run_shock(semi, LinRange(0, 0.005, 1024), 0.8, limiter!, callbacks=())

# ╔═╡ c9ae297a-5531-441d-b5d9-c5171b00b0e2
function xtdata(sol, semi; nvisnodes=nothing)
	pd = [PlotData1D(sol.u[i], semi; nvisnodes).data for i in eachindex(sol.u)]
	data = Dict(Symbol(name) => reduce(hcat, getindex.(pd, :, i))
		 for (i, name) in enumerate(varnames(cons2prim, semi.equations)))
	x = PlotData1D(sol.u[1], semi; nvisnodes).x
	t = sol.t
	(; x, t, data)
end

# ╔═╡ 44c6f87b-1f80-40ea-ba55-efac18a35aee
x, t, data = xtdata(sol, semi)

# ╔═╡ 3ce86cb1-9613-4d3c-b8f7-88138e57b070
iₓ_exit = findfirst(>(sum(last.(slabs)[1:2])), x)

# ╔═╡ 538a6c74-a16c-4241-9102-f480dbd411b2
heatmap(x, t, data[:p])

# ╔═╡ c45d8fe2-c37d-4322-8ce1-f0a2f027f03a
dt = t[2] - t[1]

# ╔═╡ 5738d94c-ea9d-4381-8985-6ed62498b89d
plot(t, data[:v1][iₓ_exit, :])

# ╔═╡ e005a639-2a04-40c6-9976-ae7757cd20cc
plot(t, cumsum(data[:v1][iₓ_exit, :] .* dt) ./ ustrip(u"m", 0.875u"inch"))

# ╔═╡ Cell order:
# ╠═f1b8cf12-7bd3-11ed-0e4d-6934b16fe50f
# ╠═b9cc394e-5951-40ae-b53d-f9454fe20421
# ╠═f3c67e3f-6b14-4493-a9e1-71754375c30b
# ╠═3cee7d89-132e-4b15-845e-e9ec510539e3
# ╠═68b78569-9c4a-4b79-b6cd-089aa53e0de1
# ╠═edc0fb1c-b9e4-4274-bf08-aba90188372b
# ╠═756d3adb-b79d-4737-9e18-0edf08b9ad5f
# ╠═3ce86cb1-9613-4d3c-b8f7-88138e57b070
# ╠═b8d1b888-d501-47a8-a0b3-670ea10ad6b0
# ╠═538a6c74-a16c-4241-9102-f480dbd411b2
# ╠═c45d8fe2-c37d-4322-8ce1-f0a2f027f03a
# ╠═5738d94c-ea9d-4381-8985-6ed62498b89d
# ╠═e005a639-2a04-40c6-9976-ae7757cd20cc
# ╠═d87dec7b-7f0e-402e-900b-a9263c894a87
# ╠═10bda30d-a84f-41d9-8c21-585f976cbcb7
# ╠═14b66293-209b-4e61-b7dc-f52755bf3378
# ╠═d0126bd1-24c3-45f4-97e7-39a4fc4e00ae
# ╠═6b0c4c1d-f81e-4357-ba24-679c2c06ded9
# ╠═95ff2fa9-d64b-4aa2-bee2-e708c64edd70
# ╠═c9ae297a-5531-441d-b5d9-c5171b00b0e2
# ╠═44c6f87b-1f80-40ea-ba55-efac18a35aee
