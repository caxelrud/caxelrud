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
  behavior — see `sanchez_lacombe.jl` — plus binary *mixture* PVT behavior
  via van der Waals-type mixing rules — see `sanchez_lacombe_mixture.jl` —
  plus mixture chemical potentials/activities, restricted to a
  thermodynamically consistent (constant local hole volume) formulation
  following von Konigslow, Park & Thompson (2017) — see
  `sanchez_lacombe_activity.jl`.
- Free-radical polymerization kinetics (QSSA rate expressions, molecular
  weight averages, isothermal batch conversion) — see
  `kinetics_free_radical.jl`.
- Step-growth (condensation) polymerization kinetics and the Flory
  molecular weight distribution — see `kinetics_step_growth.jl`.
- Coordination (Ziegler-Natta / metallocene) polymerization kinetics — see
  `kinetics_coordination.jl`.
- Ideal CSTR, PFR, and CSTR-train reactor unit operations for all three
  kinetic schemes above — see `reactors.jl`.
- A single PFR-with-recycle unit operation (one implicit recycle loop
  around one reactor), solved via closed forms and root-finding derived
  by hand for that one topology — see `reactor_recycle.jl`.
- A general flowsheet framework — [`Stream`](@ref)s, [`mix`](@ref)ers,
  [`split_stream`](@ref) splitters, and [`solve_tear`](@ref) for recycle
  convergence — that composes the reactor unit operations above into
  *arbitrary* topologies (multiple reactors, multiple or nested recycle
  loops, branch-and-remix networks), using the same tear-stream method
  standard process simulators use — see `flowsheet.jl`.

"""
module PolyRigorous

include("components.jl")
include("flory_huggins.jl")
include("phase_equilibrium.jl")
include("sanchez_lacombe.jl")
include("sanchez_lacombe_mixture.jl")
include("sanchez_lacombe_activity.jl")
include("kinetics_free_radical.jl")
include("kinetics_step_growth.jl")
include("kinetics_coordination.jl")
include("reactors.jl")
include("reactor_recycle.jl")
include("flowsheet.jl")

export Species, species, SPECIES_DB, degree_of_polymerization, segment_number

export chi_from_solubility, ln_activity_solvent, activity_solvent,
       ln_activity_polymer, activity_polymer, flory_huggins_activities,
       osmotic_pressure, chi_spinodal, critical_point

export binodal_pair, binodal_curve, spinodal_curve

export sl_eos_residual, reduced_density, density, specific_volume,
       thermal_expansion_coefficient, isothermal_compressibility

export sl_mixing_rules, sl_mixture_species

export sl_consistent_mixing_rules, sl_mixture_ln_activities, sl_mixture_activities

export initiation_rate, radical_concentration, propagation_rate,
       kinetic_chain_length, Xn_combination, Xn_disproportionation,
       Xn_mixed, monomer_concentration, conversion

export extent_reaction_external_catalyst, extent_reaction_self_catalyzed,
       carothers_Xn, flory_mole_fraction, flory_weight_fraction, flory_Xw,
       flory_PDI

export cstr_free_radical, pfr_free_radical, cstr_step_growth_external_catalyst,
       cstr_step_growth_self_catalyzed, pfr_step_growth_external_catalyst,
       pfr_step_growth_self_catalyzed

export cstr_train_free_radical, cstr_train_step_growth_external_catalyst,
       cstr_train_step_growth_self_catalyzed

export coordination_propagation_rate, transfer_rate_constant, Xn_coordination,
       monomer_concentration_coordination, conversion_coordination

export cstr_coordination, pfr_coordination, cstr_train_coordination

export pfr_free_radical_recycle, pfr_coordination_recycle,
       pfr_step_growth_external_catalyst_recycle,
       pfr_step_growth_self_catalyzed_recycle

export Stream, mix, split_stream, solve_tear

end # module PolyRigorous
