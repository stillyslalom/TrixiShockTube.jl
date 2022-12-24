module TrixiShockTube

using Trixi
using Trixi: varnames
using PyThermo
using PyThermo.ShockTube
using Unitful
using Roots
using OrdinaryDiffEq

export build_shockube_ic, calculate_Mach
export build_limiter, build_semidiscretization, build_ODE
export run_shock, xtdata

## Utilities

"""
    calculate_Mach(driver, driven, P_driver = ustrip(u"Pa", PyThermo.pressure(driver)))

Calculate the Mach number generated in a lossless 1D shock tube given the driver and driven gas 
compositions and thermodynamic states. The gas compositions and states are specified as either 
`Species` or `Mixture` objects from `PyThermo.jl`.

# Arguments
- `driver::PyThermo.Chemical`: driver gas composition
- `driven::PyThermo.Chemical`: driven gas composition

# Optional arguments
- `M₀::Float64`: initial guess for Mach number; default is 1.1

# Keyword arguments
- `P_driver::Float64`: driver gas pressure in Pa; taken from `driver` argument pressure if not specified.
"""
function calculate_Mach(driver, driven, M₀ = 1.1;
                        P_driver = ustrip(u"Pa", PyThermo.pressure(driver)))
    M = fzero(M₀) do M
         ustrip(u"Pa", PyThermo.pressure(driverpressure(driver, driven, M))) - P_driver
    end
end


"""
    locator(lengths, origin)

Create a function that returns the index of the slab containing a given position.

# Arguments
- `lengths::AbstractVector{<:Real}`: length of each slab [m]
- `origin::Real`: coordinate system origin of the shock tube in meters
     relative to the start of the first slab
"""
function locator(lengths, origin)
	locations = cumsum(lengths) .- origin
	function loc(x)
		idx = findfirst(l -> x < l, locations)
		isnothing(idx) ? length(locations) : idx
	end
end

"""
    build_shockube_ic(slabs::Tuple{<:Chemical, <:Real}...)

Create an initial condition for a shock tube problem given the gas composition and length of each slab.
The gas composition is specified as either a `Species` or `Mixture` object from `PyThermo.jl`.

# Arguments
- `slabs::Tuple{<:Chemical, <:Real}...`: gas composition and length [m] of each slab 

# Keyword arguments
- `origin::Float64`: coordinate system origin of the shock tube in meters
     relative to the start of the first slab; default is 0.0.
- `pseudodensity::Float64`: artificial fill fraction relative to true density [kg/m^3] in regions
    outside of the specified slab. Increase if numerical instabilities are encountered. Default is 1e-4.
"""
function build_shocktube_ic(slabs; origin = 0.0, pseudodensity = 1e-4)
	N = length(slabs)
	gases = first.(slabs)
	lengths = last.(slabs)
	gammas = isentropic_exponent.(gases)
	gas_constants = ustrip.(u"kJ/(kg*K)", R_specific.(gases))
	
	equations = CompressibleEulerMulticomponentEquations1D(;
		gammas, gas_constants)
	
	loc = locator(lengths, origin)
	ρs = @. ustrip(u"kg/m^3", PyThermo.density(gases))
	ps = @. ustrip(u"Pa", PyThermo.pressure(gases))
	
	ic = function initial_condition_shock(x::SVector{1,T}, t, equations) where {T}
		idx = loc(only(x))
		p = ps[idx]
		ρ = SVector{N,T}(ρs[idx]*(i == idx ? 1.0 : pseudodensity) for i in 1:N)
		return prim2cons(vcat(zero(T), p, ρ), equations)
	end
	return ic, equations
end

## Trixi.jl missing methods

"""
    boundary_condition_reflect(u_inner, orientation, direction, x, t,
                               surface_flux_function, equations::CompressibleEulerMulticomponentEquations1D)

Boundary condition reflecting the flow at the domain edges.
"""
function boundary_condition_reflect(u_inner, orientation, direction, x, t,
    surface_flux_function, equations::CompressibleEulerMulticomponentEquations1D)

  u_inner_reflect = SVector(-u_inner[1], u_inner[2:end]...)
  # Calculate boundary flux
  if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
    flux = surface_flux_function(u_inner, u_inner_reflect, orientation, equations)
  else # u_boundary is "left" of boundary, u_inner is "right" of boundary
    flux = surface_flux_function(u_inner_reflect, u_inner, orientation, equations)
  end

  return flux
end

"""
    build_limiter(thresh, equations::CompressibleEulerMulticomponentEquations1D)
    build_limiter(equations::CompressibleEulerMulticomponentEquations1D)

Build a `PositivityPreservingLimiterZhangShu` for the `CompressibleEulerMulticomponentEquations1D`
equations with the given threshold(s) for each of N components' (`NC`) density and the total pressure.

# Arguments
- `thresh::NTuple{NC+1, <:Real}`: minimum threshold(s) for each species' density and the total pressure
- `equations::CompressibleEulerMulticomponentEquations1D`: equations for which to build the limiter

Alternatively, a single threshold can be specified for all species' densities and the total pressure.
If a threshold is not specified, a default value of 1e-7 is used.
"""
function build_limiter(thresh, equations::CompressibleEulerMulticomponentEquations1D{NV,NC,T}) where {NV,NC,T}
    variables = (ntuple(i -> (u, e) -> u[i+2], NC)..., Trixi.pressure)
    thresholds = ntuple(i -> T(thresh[i]), NC+1)
    PositivityPreservingLimiterZhangShu(; thresholds, variables)
end

function build_limiter(thresh::Real, equations::CompressibleEulerMulticomponentEquations1D{NV,NC,T}) where {NV,NC,T}
    thresholds = ntuple(i -> thresh, NC+1)
    build_limiter(thresholds, equations)
end

function build_limiter(equations::CompressibleEulerMulticomponentEquations1D{NV,NC,T}) where {NV,NC,T}
    build_limiter(1e-7, equations)
end

"""
    density_pressure(u, equations::CompressibleEulerMulticomponentEquations1D)

Compute the density times the pressure from the conserved variables `u`.
"""
@inline function density_pressure(u, equations::CompressibleEulerMulticomponentEquations1D)
    rho_v1, rho_e = u
  
    rho          = Trixi.density(u, equations)
    gamma        = Trixi.totalgamma(u, equations)
    rho_times_p  = (gamma - 1) * (rho * rho_e - 0.5 * (rho_v1^2))
  
    return rho_times_p
  end

## Initialization

"""
    build_semidiscretization(slabs, equations::CompressibleEulerMulticomponentEquations1D;
                             origin = 0.0,
                             initial_refinement_level = 7,
                             n_cells_max = 10_000,
                             surface_flux = flux_lax_friedrichs,
                             volume_flux  = flux_ranocha,
                             basis        = LobattoLegendreBasis(3),
                             indicator_sc = IndicatorHennemannGassner(equations, basis,
                                                                     alpha_max=0.5,
                                                                     alpha_min=0.001,
                                                                     alpha_smooth=true,
                                                                     variable=density_pressure))

Build a `DGSEM` semidiscretization with the specified `surface_flux` and `volume_flux` using the
`CompressibleEulerMulticomponentEquations1D` equations.

# Arguments
- `slabs::Tuple{<:Chemical, <:Real}...`: gas composition and length [m] of each slab
- `equations::CompressibleEulerMulticomponentEquations1D`: equations for which to build the semidiscretization

# Keyword arguments
- `origin::Real`: origin of the domain
- `initial_refinement_level::Integer`: initial refinement level (number of elements = 2^level)
- `n_cells_max::Integer`: maximum number of cells
- `surface_flux::Function`: surface flux function
- `volume_flux::Function`: two-point volume flux function (must be one of
    `flux_ranocha`, `flux_shima_etal`, `flux_chandrashekar` or `flux_kennedy_gruber`)
- `basis::LobattoLegendreBasis`: polynomial basis of the approximation space of specified order
- `indicator_sc::IndicatorHennemannGassner`: shock-capturing indicator for the DGSEM solver
"""
function build_semidiscretization(slabs, ic, equations::CompressibleEulerMulticomponentEquations1D;
        origin       = 0.0,
        initial_refinement_level = 7,
        n_cells_max = 10_000,
        surface_flux = flux_lax_friedrichs,
        volume_flux  = flux_ranocha,
        basis        = LobattoLegendreBasis(3),
        indicator_sc = IndicatorHennemannGassner(equations, basis,
                                                 alpha_max=0.5,
                                                 alpha_min=0.001,
                                                 alpha_smooth=true,
                                                 variable=density_pressure))

    volume_integral     = VolumeIntegralShockCapturingHG(indicator_sc;
                                                        volume_flux_dg=volume_flux,
                                                        volume_flux_fv=surface_flux)
    solver              = DGSEM(basis, surface_flux, volume_integral)

    coordinates_min = (-origin,)
    coordinates_max = (sum(last.(slabs)) - origin,)
    mesh = TreeMesh(coordinates_min, coordinates_max;
                initial_refinement_level,
                n_cells_max,
                periodicity=false)

    semi = SemidiscretizationHyperbolic(mesh, equations, ic, solver; boundary_conditions=boundary_condition_reflect)
end


"""
    build_ODE(semi, saveat, cfl; callbacks)

Build an ODE problem with the specified `semi` semidiscretization, `saveat` output time steps and `cfl` CFL limiter
using the `CompressibleEulerMulticomponentEquations1D` equations. Any additional `callbacks` are passed to the
`CallbackSet` constructor.

# Arguments
- `semi::SemidiscretizationHyperbolic`: semidiscretization for which to build the ODE problem
- `saveat::AbstractVector{<:Real}`: time steps at which to save the solution
- `cfl::Real`: CFL limiter

# Keyword arguments
- `callbacks::Tuple(Function)...`: additional callbacks to pass to the `CallbackSet` constructor. See
    [Trixi.jl callbacks](https://trixi-framework.github.io/Trixi.jl/stable/callbacks/) for more details.
"""
function build_ODE(semi, saveat, cfl; callbacks=())
    tspan = extrema(saveat)

	ode = semidiscretize(semi, tspan)

	stepsize_callback = StepsizeCallback(; cfl)

	callback = CallbackSet(callbacks...,
                        stepsize_callback)

    return ode, callback
end

"""
    run_shock(semi, saveat, cfl, limiter!; callbacks)

Run a shock tube simulation with the specified `semi` semidiscretization, `saveat` output time steps, `cfl` CFL limiter,
and `limiter!` positivity-preserving limiter using the `CompressibleEulerMulticomponentEquations1D` equations. Any additional
`callbacks` are passed to the `CallbackSet` constructor.

# Arguments
- `semi::SemidiscretizationHyperbolic`: semidiscretization for which to build the ODE problem
- `saveat::AbstractVector{<:Real}`: time steps at which to save the solution
- `cfl::Real`: CFL limiter
- `limiter!::Function`: positivity-preserving limiter

# Keyword arguments
- `callbacks::Tuple(Function)...`: additional callbacks to pass to the `CallbackSet` constructor. See
    [Trixi.jl callbacks](https://trixi-framework.github.io/Trixi.jl/stable/callbacks/) for more details.
"""
function run_shock(semi, saveat, cfl, limiter!; callbacks=())
    ode, callback = build_ODE(semi, saveat, cfl; callbacks)
    solve(ode, 
    CarpenterKennedy2N54(limiter!, limiter!, williamson_condition=false);
          dt=1.0, # arbitrary value
		  saveat,
          save_everystep=false,
          callback)
end

"""
    xtdata(sol, semi; nvisnodes=nothing)

Extract the primitive field solution data from the `sol` solution and `semi` semidiscretization at `nvisnodes` visualization nodes.
By default, `nvisnodes` is set to twice the number of nodes in the polynomial approximation space. If `nvisnodes`
is set to zero, the solution data is extracted at the nodes of the approximation space.

# Arguments
- `sol::TrixiODESolution`: solution to extract the data from
- `semi::SemidiscretizationHyperbolic`: semidiscretization for which to extract the data

# Keyword arguments
- `nvisnodes::Integer`: number of visualization nodes
"""
function xtdata(sol, semi; nvisnodes=nothing)
	pd = [PlotData1D(sol.u[i], semi; nvisnodes).data for i in eachindex(sol.u)]
	data = Dict(Symbol(name) => reduce(hcat, getindex.(pd, :, i))
		 for (i, name) in enumerate(varnames(cons2prim, semi.equations)))
	x = PlotData1D(sol.u[1], semi; nvisnodes).x
	t = sol.t
	(; x, t, data)
end

end # module