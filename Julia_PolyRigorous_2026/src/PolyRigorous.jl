"""
    PolyRigorous

A Julia toolkit for rigorous polymer thermodynamics and kinetics — the kind
of calculations underlying tools like Aspen Polymers Plus — covering:

- A small built-in species database ([`species`](@ref), [`SPECIES_DB`](@ref)).
- Flory-Huggins lattice theory (activities, osmotic pressure, spinodal and
  critical point) — see `flory_huggins.jl`.
- Liquid-liquid binodal/spinodal phase-split calculations built on top of
  Flory-Huggins — see `phase_equilibrium.jl`.
- The Sanchez-Lacombe lattice-fluid equation of state for pure-component PVT
  behavior — see `sanchez_lacombe.jl`.
- Free-radical polymerization kinetics (QSSA rate expressions, molecular
  weight averages, isothermal batch conversion) — see
  `kinetics_free_radical.jl`.
- Step-growth (condensation) polymerization kinetics and the Flory
  molecular weight distribution — see `kinetics_step_growth.jl`.
- Ideal CSTR and PFR reactor unit operations built on both kinetic
  schemes — see `reactors.jl`.

Full flowsheet simulation (connecting multiple unit operations with
streams and recycle) is still out of scope for this version. See the
repository README for the roadmap.
"""
module PolyRigorous

include("components.jl")
include("flory_huggins.jl")
include("phase_equilibrium.jl")
include("sanchez_lacombe.jl")
include("kinetics_free_radical.jl")
include("kinetics_step_growth.jl")
include("reactors.jl")

export Species, species, SPECIES_DB, degree_of_polymerization, segment_number

export chi_from_solubility, ln_activity_solvent, activity_solvent,
       ln_activity_polymer, activity_polymer, flory_huggins_activities,
       osmotic_pressure, chi_spinodal, critical_point

export binodal_pair, binodal_curve, spinodal_curve

export sl_eos_residual, reduced_density, density, specific_volume,
       thermal_expansion_coefficient, isothermal_compressibility

export initiation_rate, radical_concentration, propagation_rate,
       kinetic_chain_length, Xn_combination, Xn_disproportionation,
       Xn_mixed, monomer_concentration, conversion

export extent_reaction_external_catalyst, extent_reaction_self_catalyzed,
       carothers_Xn, flory_mole_fraction, flory_weight_fraction, flory_Xw,
       flory_PDI

export cstr_free_radical, pfr_free_radical, cstr_step_growth_external_catalyst,
       cstr_step_growth_self_catalyzed, pfr_step_growth_external_catalyst,
       pfr_step_growth_self_catalyzed

end # module PolyRigorous
