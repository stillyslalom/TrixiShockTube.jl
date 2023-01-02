var documenterSearchIndex = {"docs":
[{"location":"reference/","page":"Function reference","title":"Function reference","text":"","category":"page"},{"location":"reference/","page":"Function reference","title":"Function reference","text":"Modules = [TrixiShockTube]","category":"page"},{"location":"reference/#TrixiShockTube.boundary_condition_reflect-Tuple{Any, Any, Any, Any, Any, Any, Trixi.CompressibleEulerMulticomponentEquations1D}","page":"Function reference","title":"TrixiShockTube.boundary_condition_reflect","text":"boundary_condition_reflect(u_inner, orientation, direction, x, t,\n                           surface_flux_function, equations::CompressibleEulerMulticomponentEquations1D)\n\nBoundary condition reflecting the flow at the domain edges.\n\n\n\n\n\n","category":"method"},{"location":"reference/#TrixiShockTube.build_ODE-Tuple{Any, Any, Any}","page":"Function reference","title":"TrixiShockTube.build_ODE","text":"build_ODE(semi, saveat, cfl; callbacks)\n\nBuild an ODE problem with the specified semi semidiscretization, saveat output time steps and cfl CFL limiter using the CompressibleEulerMulticomponentEquations1D equations. Any additional callbacks are passed to the CallbackSet constructor.\n\nArguments\n\nsemi::SemidiscretizationHyperbolic: semidiscretization for which to build the ODE problem\nsaveat::AbstractVector{<:Real}: time steps at which to save the solution\ncfl::Real: CFL limiter\n\nKeyword arguments\n\ncallbacks::Tuple(Function)...: additional callbacks to pass to the CallbackSet constructor. See   Trixi.jl callbacks for more details.\n\n\n\n\n\n","category":"method"},{"location":"reference/#TrixiShockTube.build_limiter-Union{Tuple{T}, Tuple{NC}, Tuple{NV}, Tuple{Any, Trixi.CompressibleEulerMulticomponentEquations1D{NV, NC, T}}} where {NV, NC, T}","page":"Function reference","title":"TrixiShockTube.build_limiter","text":"build_limiter(thresh, equations::CompressibleEulerMulticomponentEquations1D)\nbuild_limiter(equations::CompressibleEulerMulticomponentEquations1D)\n\nBuild a PositivityPreservingLimiterZhangShu for the CompressibleEulerMulticomponentEquations1D equations with the given threshold(s) for each of N components' (NC) density and the total pressure.\n\nArguments\n\nthresh::NTuple{NC+1, <:Real}: minimum threshold(s) for each species' density and the total pressure\nequations::CompressibleEulerMulticomponentEquations1D: equations for which to build the limiter\n\nAlternatively, a single threshold can be specified for all species' densities and the total pressure. If a threshold is not specified, a default value of 1e-7 is used.\n\n\n\n\n\n","category":"method"},{"location":"reference/#TrixiShockTube.build_semidiscretization-Tuple{Any, Any, Trixi.CompressibleEulerMulticomponentEquations1D}","page":"Function reference","title":"TrixiShockTube.build_semidiscretization","text":"build_semidiscretization(slabs, equations::CompressibleEulerMulticomponentEquations1D;\n                         origin = 0.0,\n                         initial_refinement_level = 7,\n                         n_cells_max = 10_000,\n                         surface_flux = flux_lax_friedrichs,\n                         volume_flux  = flux_ranocha,\n                         basis        = LobattoLegendreBasis(3),\n                         indicator_sc = IndicatorHennemannGassner(equations, basis,\n                                                                 alpha_max=0.5,\n                                                                 alpha_min=0.001,\n                                                                 alpha_smooth=true,\n                                                                 variable=density_pressure))\n\nBuild a DGSEM semidiscretization with the specified surface_flux and volume_flux using the CompressibleEulerMulticomponentEquations1D equations.\n\nArguments\n\nslabs::Tuple{<:Chemical, <:Real}...: gas composition and length [m] of each slab\nequations::CompressibleEulerMulticomponentEquations1D: equations for which to build the semidiscretization\n\nKeyword arguments\n\norigin::Real: origin of the domain\ninitial_refinement_level::Integer: initial refinement level (number of elements = 2^level)\nn_cells_max::Integer: maximum number of cells\nsurface_flux::Function: surface flux function\nvolume_flux::Function: two-point volume flux function (must be one of   flux_ranocha, flux_shima_etal, flux_chandrashekar or flux_kennedy_gruber)\nbasis::LobattoLegendreBasis: polynomial basis of the approximation space of specified order\nindicator_sc::IndicatorHennemannGassner: shock-capturing indicator for the DGSEM solver\n\n\n\n\n\n","category":"method"},{"location":"reference/#TrixiShockTube.build_shocktube_ic-Tuple{Any}","page":"Function reference","title":"TrixiShockTube.build_shocktube_ic","text":"build_shockube_ic(slabs::Tuple{<:Chemical, <:Real}...)\n\nCreate an initial condition for a shock tube problem given the gas composition and length of each slab. The gas composition is specified as either a Species or Mixture object from PyThermo.jl.\n\nArguments\n\nslabs::Tuple{<:Chemical, <:Real}...: gas composition and length [m] of each slab \n\nKeyword arguments\n\norigin::Float64: coordinate system origin of the shock tube in meters    relative to the start of the first slab; default is 0.0.\npseudodensity::Float64: artificial fill fraction relative to true density [kg/m^3] in regions   outside of the specified slab. Increase if numerical instabilities are encountered. Default is 1e-4.\n\n\n\n\n\n","category":"method"},{"location":"reference/#TrixiShockTube.calculate_Mach","page":"Function reference","title":"TrixiShockTube.calculate_Mach","text":"calculate_Mach(driver, driven, P_driver = ustrip(u\"Pa\", PyThermo.pressure(driver)))\n\nCalculate the Mach number generated in a lossless 1D shock tube given the driver and driven gas  compositions and thermodynamic states. The gas compositions and states are specified as either  Species or Mixture objects from PyThermo.jl.\n\nArguments\n\ndriver::PyThermo.Chemical: driver gas composition\ndriven::PyThermo.Chemical: driven gas composition\n\nOptional arguments\n\nM₀::Float64: initial guess for Mach number; default is 1.1\n\nKeyword arguments\n\nP_driver::Float64: driver gas pressure in Pa; taken from driver argument pressure if not specified.\n\n\n\n\n\n","category":"function"},{"location":"reference/#TrixiShockTube.density_pressure-Tuple{Any, Trixi.CompressibleEulerMulticomponentEquations1D}","page":"Function reference","title":"TrixiShockTube.density_pressure","text":"density_pressure(u, equations::CompressibleEulerMulticomponentEquations1D)\n\nCompute the density times the pressure from the conserved variables u.\n\n\n\n\n\n","category":"method"},{"location":"reference/#TrixiShockTube.locator-Tuple{Any, Any}","page":"Function reference","title":"TrixiShockTube.locator","text":"locator(lengths, origin)\n\nCreate a function that returns the index of the slab containing a given position.\n\nArguments\n\nlengths::AbstractVector{<:Real}: length of each slab [m]\norigin::Real: coordinate system origin of the shock tube in meters    relative to the start of the first slab\n\n\n\n\n\n","category":"method"},{"location":"reference/#TrixiShockTube.run_shock-NTuple{4, Any}","page":"Function reference","title":"TrixiShockTube.run_shock","text":"run_shock(semi, saveat, cfl, limiter!; callbacks)\n\nRun a shock tube simulation with the specified semi semidiscretization, saveat output time steps, cfl CFL limiter, and limiter! positivity-preserving limiter using the CompressibleEulerMulticomponentEquations1D equations. Any additional callbacks are passed to the CallbackSet constructor.\n\nArguments\n\nsemi::SemidiscretizationHyperbolic: semidiscretization for which to build the ODE problem\nsaveat::AbstractVector{<:Real}: time steps at which to save the solution\ncfl::Real: CFL limiter\nlimiter!::Function: positivity-preserving limiter\n\nKeyword arguments\n\ncallbacks::Tuple(Function)...: additional callbacks to pass to the CallbackSet constructor. See   Trixi.jl callbacks for more details.\n\n\n\n\n\n","category":"method"},{"location":"reference/#TrixiShockTube.xtdata-Tuple{Any, Any}","page":"Function reference","title":"TrixiShockTube.xtdata","text":"xtdata(sol, semi; nvisnodes=nothing)\n\nExtract the primitive field solution data from the sol solution and semi semidiscretization at nvisnodes visualization nodes. By default, nvisnodes is set to twice the number of nodes in the polynomial approximation space. If nvisnodes is set to zero, the solution data is extracted at the nodes of the approximation space.\n\nArguments\n\nsol::TrixiODESolution: solution to extract the data from\nsemi::SemidiscretizationHyperbolic: semidiscretization for which to extract the data\n\nKeyword arguments\n\nnvisnodes::Integer: number of visualization nodes\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = TrixiShockTube","category":"page"},{"location":"#TrixiShockTube","page":"Home","title":"TrixiShockTube","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"TrixiShockTube provides utility functions to facilitate simulations of multicomponent 1D shock tube problems using Trixi.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Let's start with a simple shock tube problem consisting of nitrogen in the driver section pressurized to 3 bar, and a 50/50 mixture of nitrogen and argon at standard temperature and pressure in the driven section. Gas properties and dimensional units are provided by PyThermo and Unitful, respectively.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using TrixiShockTube, Unitful, PyThermo\nN₂ = Species(\"N2\"; P = 3.0u\"bar\")\nN₂Ar = Mixture([\"N2\" => 0.5, \"Ar\" => 0.5])\ndensity(N₂), density(N₂Ar)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The shock tube is initialized at rest with a discontinuity at the interface between the driver and driven sections.  build_shocktube_ic creates an initializer and defines the constitutive equations from a list of Species/Mixtures and  the lengths of the corresponding shock tube sections. The shock tube is then discretized using a uniform mesh with 2^8 elements.","category":"page"},{"location":"","page":"Home","title":"Home","text":"slabs = (N₂   => ustrip(u\"m\", 30.0u\"inch\"),\n\t\tN₂Ar => ustrip(u\"m\", 125.0u\"inch\"))\nic, equations = TrixiShockTube.build_shocktube_ic(slabs)\nsemi = build_semidiscretization(slabs, ic, equations;\n\tinitial_refinement_level=8)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The shock tube is then evolved in time using a 5-stage, 4th order Runge-Kutta method restricted by a CFL number of 0.5. A positivity-preserving limiter is applied to the solution at each time step. The solution is saved at 1000 equally spaced time steps. Trixi.jl saves the conserved variables, so we need to convert the solution to primitive variables u1, p, rho1, rho2, ... with xtdata before plotting.","category":"page"},{"location":"","page":"Home","title":"Home","text":"saveat = LinRange(0, 0.02, 1000)\nlimiter! = build_limiter(equations)\nsol = run_shock(semi, saveat, 0.5, limiter!);\nx, t, data = xtdata(sol, semi)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The solution is then plotted using CairoMakie","category":"page"},{"location":"","page":"Home","title":"Home","text":"using CairoMakie\nfig, ax, hm = heatmap(x, t, data[:p];\n    axis = (xlabel = \"x [m]\", ylabel = \"t [s]\"));\ncb = Colorbar(fig[:, end+1], hm; label = \"Pressure [Pa]\");\nsave(\"pressure.png\", fig);","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: )","category":"page"},{"location":"","page":"Home","title":"Home","text":"Let's compare to the analytic solution for the density in the reflected shock region.  The analytic solution is provided by PyThermo's unexported ShockTube module. The Mach  number of the incident shock is calculated using calculate_Mach from TrixiShockTube.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using PyThermo.ShockTube\nMₛ = calculate_Mach(N₂, N₂Ar)\nsc = shockcalc(N₂, N₂Ar, Mₛ)\ndensity(sc.reflected) - maximum(data[:rho2][end, :])*u\"kg/m^3\"","category":"page"}]
}
