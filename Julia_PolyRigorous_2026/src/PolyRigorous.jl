"""
    PolyRigorous

A Julia toolkit for rigorous polymer solution thermodynamics — the kind of
calculations underlying tools like Aspen Polymers Plus — starting with:

- A small built-in species database ([`species`](@ref), [`SPECIES_DB`](@ref)).
- Flory-Huggins lattice theory (activities, osmotic pressure, spinodal and
  critical point) — see `flory_huggins.jl`.
- Liquid-liquid binodal/spinodal phase-split calculations built on top of
  Flory-Huggins — see `phase_equilibrium.jl`.
- The Sanchez-Lacombe lattice-fluid equation of state for pure-component PVT
  behavior — see `sanchez_lacombe.jl`.

This is deliberately a *thermodynamics core*: polymerization kinetics,
reactor models, and full flowsheet simulation are out of scope for this
version. See the repository README for the roadmap.
"""
module PolyRigorous

include("components.jl")
include("flory_huggins.jl")
include("phase_equilibrium.jl")
include("sanchez_lacombe.jl")

export Species, species, SPECIES_DB, degree_of_polymerization, segment_number

export chi_from_solubility, ln_activity_solvent, activity_solvent,
       ln_activity_polymer, activity_polymer, flory_huggins_activities,
       osmotic_pressure, chi_spinodal, critical_point

export binodal_pair, binodal_curve, spinodal_curve

export sl_eos_residual, reduced_density, density, specific_volume,
       thermal_expansion_coefficient, isothermal_compressibility

end # module PolyRigorous
