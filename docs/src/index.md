```@meta
CurrentModule = TrixiShockTube
```

# TrixiShockTube

[TrixiShockTube](https://github.com/stillyslalom/TrixiShockTube.jl) provides utility functions
to facilitate simulations of multicomponent 1D shock tube problems using [Trixi](https://github.com/trixi-framework/Trixi.jl).

Let's start with a simple shock tube problem consisting of nitrogen in the driver section pressurized to 3 bar, and a 50/50 mixture of nitrogen and argon at standard temperature and pressure in the driven section. Gas properties and dimensional units are provided by
[PyThermo](https://github.com/stillyslalom/PyThermo.jl) and [Unitful](https://github.com/PainterQubits/Unitful.jl), respectively.
```@repl 1
using TrixiShockTube, Unitful, PyThermo
N₂ = Species("N2"; P = 3.0u"bar")
N₂Ar = Mixture(["N2" => 0.5, "Ar" => 0.5])
density(N₂), density(N₂Ar)
```

The shock tube is initialized at rest with a discontinuity at the interface between the driver and driven sections. 
`build_shocktube_ic` creates an initializer and defines the constitutive equations from a list of `Species`/`Mixture`s and 
the lengths of the corresponding shock tube sections. The shock tube is then discretized using a uniform mesh with 2^8 elements.
```@repl 1
slabs = (N₂   => ustrip(u"m", 30.0u"inch"),
		 N₂Ar => ustrip(u"m", 125.0u"inch"))
ic, equations = TrixiShockTube.build_shocktube_ic(slabs)
semi = build_semidiscretization(slabs, ic, equations;
	initial_refinement_level=8)
```

The shock tube is then evolved in time using a 5-stage, 4th order Runge-Kutta method restricted by a CFL number of 0.5.
A positivity-preserving limiter is applied to the solution at each time step. The solution is saved at 1000 equally spaced time steps.
Trixi.jl saves the primitive variables, so we need to convert the solution to conserved variables with `xtdata` before plotting.
```@repl 1
saveat = LinRange(0, 0.02, 1000)
limiter! = build_limiter(equations)
sol = run_shock(semi, saveat, 0.5, limiter!);
x, t, data = xtdata(sol, semi)
```

The solution is then plotted using [CairoMakie](https://github.com/MakieOrg/Makie.jl/tree/master/CairoMakie)
```@repl 1
using CairoMakie
fig, ax, hm = heatmap(x, t, data[:p];
    axis = (xlabel = "x [m]", ylabel = "t [s]"));
cb = Colorbar(fig[:, end+1], hm; label = "Pressure [Pa]");
save("pressure.png", fig);
```
![](pressure.png)

Let's compare to the analytic solution for the density in the reflected shock region. 
The analytic solution is provided by PyThermo's unexported `ShockTube` module. The Mach 
number of the incident shock is calculated using `calculate_Mach` from TrixiShockTube.
```@repl 1
using PyThermo.ShockTube
Mₛ = calculate_Mach(N₂, N₂Ar)
sc = shockcalc(N₂, N₂Ar, Mₛ)
density(sc.reflected) - maximum(data[:rho2][end, :])*u"kg/m^3"
```
