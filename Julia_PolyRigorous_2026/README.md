# PolyRigorous

A Julia + [Pluto.jl](https://plutojl.org/) toolkit for rigorous polymer
solution thermodynamics — the kind of calculations underlying commercial
tools like Aspen Polymers Plus. It's a from-scratch, open implementation,
not a clone of any proprietary source or model — see [Scope & honesty
about the data](#scope--honesty-about-the-data) below.

## What's here (v0.1.0)

This first version covers **thermodynamics and kinetics**:

- A small built-in species database (`species`, `SPECIES_DB`) with a
  handful of common polymers — polystyrene, PMMA, polyethylene (both an
  HDPE-like grade and, separately, LDPE with parameters regressed from
  real PVT data in von Konigslow, Thompson et al., *Soft Matter* **14**,
  4603 (2018)), and polypropylene (isotactic, linear, same source) — and
  solvents (toluene, benzene, cyclohexane, n-hexane).
- **Flory-Huggins lattice theory**: solvent/polymer activities, osmotic
  pressure, the spinodal curve, and the critical point, for a binary
  polymer-solvent system (`src/flory_huggins.jl`).
- **Liquid-liquid phase equilibrium**: binodal (coexistence) curve
  computation by solving the equal-activity conditions between two phases
  (`src/phase_equilibrium.jl`).
- **Sanchez-Lacombe lattice-fluid equation of state** for pure-component
  PVT behavior: density, specific volume, thermal expansion coefficient,
  isothermal compressibility (`src/sanchez_lacombe.jl`); plus binary
  **mixture** PVT behavior (mixture density) via van der Waals-type
  mixing rules for the characteristic parameters, reusing the same
  pure-component EOS solver (`src/sanchez_lacombe_mixture.jl`); plus
  binary **mixture chemical potentials/activities**
  (`sl_mixture_activities`), restricted to a thermodynamically consistent
  formulation (constant local hole volume, following von Konigslow, Park &
  Thompson, *Phys. Rev. Applied* **8**, 044009 (2017)) with its own,
  separate mixing rules verified against the Euler relation
  (`src/sanchez_lacombe_activity.jl`) — see [Scope & honesty about the
  data](#scope--honesty-about-the-data) for why this needed a different
  mixing-rule convention from the PVT module above.
- **Free-radical polymerization kinetics**: QSSA rate expressions,
  kinetic chain length, number-average degree of polymerization
  (combination / disproportionation / mixed termination), and isothermal
  batch conversion vs. time (`src/kinetics_free_radical.jl`).
- **Step-growth (condensation) polymerization kinetics**: extent of
  reaction for externally-catalyzed and self-catalyzed kinetics, the
  Carothers equation (including stoichiometric imbalance), and the Flory
  "most probable" molecular weight distribution
  (`src/kinetics_step_growth.jl`).
- **Coordination (Ziegler-Natta / metallocene) polymerization kinetics**:
  the "living-like" chain-growth mechanism, where Xₙ is derived from
  conservation of chain count (propagation rate over the total rate of
  chain-releasing events) rather than an uncertain closed-form formula,
  plus the classic industrial hydrogen-response lever for molecular
  weight control (`src/kinetics_coordination.jl`).
- **Reactor unit operations** built on all three kinetic schemes above:
  ideal CSTR (steady-state mass balance), PFR (kinetically equivalent to a
  batch reactor run for a time equal to its residence time), and
  CSTRs-in-series reactor trains — verified against the classic result
  that a train converges to PFR performance as the stage count grows at
  fixed total residence time (`src/reactors.jl`).
- **PFR with recycle**: a single ideal PFR whose outlet is split, with a
  fraction recycled back and mixed into the fresh feed, for all three
  kinetic schemes. The recycle stream's composition depends on the very
  outlet it feeds, making this a genuinely implicit ("flowsheet-style")
  problem rather than a single forward calculation — solved in closed
  form for the free-radical and coordination cases and the step-growth
  external-catalyst case, and by bracketed root-finding for the
  step-growth self-catalyzed case. Reduces exactly to the plain PFR at
  recycle ratio `R = 0`, and approaches the corresponding CSTR (same
  nominal residence time) as `R → ∞` — the classic "PFR with infinite
  recycle behaves like a CSTR" result, checked directly in the test suite
  (`src/reactor_recycle.jl`).
- **A general flowsheet framework**: `Stream`s carrying a volumetric flow
  rate and composition, `mix` and `split_stream` unit operations, and the
  reactor functions above extended (via new methods, same names) to
  consume and produce `Stream`s — so arbitrary topologies (multiple
  reactors, branch-and-remix networks, multiple or nested recycle loops)
  are just ordinary Julia function composition. Recycle loops close by
  *tear-stream* convergence (`solve_tear`): guess the torn stream, go
  around the loop once, iterate the guess toward its own image by damped
  successive substitution — the standard "sequential-modular" method used
  by commercial process simulators, not something specific to this
  package's own derivations. Verified by construction: assembling the
  *same* single-PFR-with-recycle topology this way and solving it with
  `solve_tear` reproduces every one of `reactor_recycle.jl`'s hand-derived
  closed forms to numerical precision (`src/flowsheet.jl`).
- **Multi-zone reactors with staged initiator injection**: `n` zones/segments
  in series, each with *additional* fresh initiator injected on top of
  whatever survives (decays) from the previous stage — the standard model
  for the two high-pressure (150-300 MPa) LDPE free-radical process
  configurations, a multi-zone **autoclave**
  (`cstr_train_free_radical_staged`, each zone an ideal CSTR) and a
  **tubular** reactor with multiple injection points
  (`pfr_train_free_radical_staged`, each segment a PFR). Reduces exactly
  to a plain (non-staged) CSTR train/PFR when all the initiator is charged
  in the first stage and nothing thereafter — checked in the test suite,
  along with the property that splitting a plain PFR into extra
  zero-injection segments changes nothing (isolating the staging
  machinery from the injection feature). Also includes
  `ideal_gas_concentration`, converting a monomer partial pressure to the
  concentration this package's kinetics expect, for modeling gas-phase
  reactors (`src/reactor_staged_injection.jl`).
- An interactive **Pluto notebook**, `notebooks/ThermoExplorer.jl`, that
  puts sliders and dropdowns on top of all of the above.

## Industrial polyolefin processes → which model to use

Polyethylene and polypropylene are made industrially by several distinct
reactor/process configurations. None of them need new kinetics beyond
what's already in `kinetics_free_radical.jl`/`kinetics_coordination.jl` —
each maps onto an existing (or, for the two LDPE cases, newly added)
reactor model:

| Process | Kinetics | Model | Typical conditions (illustrative) |
| --- | --- | --- | --- |
| Gas-phase fluidized bed (e.g. UNIPOL-type PE/PP) | Coordination | `cstr_coordination` (well-mixed bed) — feed monomer via `ideal_gas_concentration(P, T)` | ~80-100 °C, monomer partial pressure ~1-2 MPa, τ ~1-4 h |
| Slurry loop (e.g. Phillips-type HDPE) | Coordination | `cstr_coordination` (the standard well-mixed idealization for a high-recirculation loop) or, more rigorously, `pfr_coordination_recycle` at large `R` (shown to converge to the same CSTR result) | ~85-110 °C, ~3-4 MPa, τ ~1-1.5 h |
| Solution (e.g. Dow-type PE/EP elastomers) | Coordination | `cstr_coordination`/`cstr_train_coordination` at high `T` | ~150-250 °C, ~4-15 MPa, τ ~1-5 min (high catalyst activity) |
| High-pressure autoclave (LDPE) | Free-radical | `cstr_train_free_radical_staged` (multi-zone, staged initiator) | ~150-300 MPa, ~150-300 °C, τ per zone ~30-60 s |
| High-pressure tubular (LDPE) | Free-radical | `pfr_train_free_radical_staged` (staged initiator injection points) | ~150-300 MPa, ~150-300 °C, τ ~1-3 min total |

The conditions above are illustrative orders of magnitude for context,
not regression targets — see [Scope & honesty about the
data](#scope--honesty-about-the-data). Section 14 of the notebook works
through gas-phase/slurry-loop/solution with the built-in polyethylene and
polypropylene species; section 15 works through the LDPE autoclave and
tubular models.

## Roadmap (not yet implemented)

Nothing from the original roadmap remains: recycle loops, coordination
reactors, a general flowsheet solver, and (most recently) Sanchez-Lacombe
mixture chemical potentials/activities are all implemented above. Future
directions would be new scope beyond what this package originally set out
to cover (e.g. a general multi-component, N-ary flowsheet solver beyond
today's binary-mixture and single/multi-reactor-with-recycle scope, or
non-isothermal reactors) rather than unfinished items.

## Getting started

You need a working [Julia](https://julialang.org/) installation
(1.9+recommended).

```bash
git clone <this repo>
cd Julia_PolyRigorous_2026
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

This has been run end-to-end (Julia 1.13, all 226 tests passing) as part of building this package.

### Using the package directly

```julia
using Pkg
Pkg.activate(".")
using PolyRigorous

poly = species("polystyrene")
solv = species("toluene")
N = degree_of_polymerization(poly, solv)

χ = chi_from_solubility(solv.Vm, solv.δ, poly.δ, 298.15)
crit = critical_point(N)
@show χ, crit
```

### Running the Pluto notebook

```julia
using Pkg
Pkg.add("Pluto")
using Pluto
Pluto.run(notebook = "notebooks/ThermoExplorer.jl")
```

The notebook activates this repository's own environment (`Pkg.activate`
pointed at the repo root) rather than Pluto's own per-notebook package
manager, so it picks up whatever you've already `Pkg.add`ed above.

## Naming convention

Julia allows Unicode identifiers, and this package uses that on purpose:
physical quantities are spelled the way they appear in the docstrings'
math, not transliterated into ASCII. So the code reads `χ` (`\chi`), `φ₂`
(`\varphi\_2`), `δ` (`\delta`), `ν` (`\nu`), and reduced Sanchez-Lacombe
variables as `T̃`, `P̃`, `ρ̃` (a letter plus a combining tilde, `\tilde`),
rather than `chi`, `phi2`, `delta`, `nu`, `T_tilde`. Rate constants follow
the textbook's own subscript convention, `k_p`, `k_d`, `k_t`. Exported
*function* names stay plain English (`chi_spinodal`, `critical_point`,
`kinetic_chain_length`, ...) so the API stays easy to type, grep, and
tab-complete without a LaTeX input method — only the parameters and local
variables *inside* those functions (and the fields of the values they
return, like `critical_point`'s `φ₂_c`/`χ_c`) use the symbols. In the
Julia REPL or an editor with LaTeX tab-completion (VS Code, Pluto, the
Julia REPL itself), typing `\chi<TAB>` produces `χ`, `\varphi<TAB>\_2<TAB>`
produces `φ₂`, and `\tilde<TAB>` after a letter adds the combining tilde.

## Project layout

```
Project.toml            # PolyRigorous package manifest
src/
  PolyRigorous.jl        # module entry point
  components.jl           # Species type + built-in database
  flory_huggins.jl         # Flory-Huggins activities, spinodal, critical point
  phase_equilibrium.jl     # binodal curve via NLsolve
  sanchez_lacombe.jl       # Sanchez-Lacombe pure-component EOS
  sanchez_lacombe_mixture.jl # Sanchez-Lacombe binary mixture PVT
  sanchez_lacombe_activity.jl # Sanchez-Lacombe mixture chemical potentials
  kinetics_free_radical.jl # free-radical polymerization kinetics (QSSA)
  kinetics_step_growth.jl  # step-growth kinetics + Flory MWD
  kinetics_coordination.jl # coordination (Ziegler-Natta) kinetics
  reactors.jl              # ideal CSTR, PFR, and CSTR-train unit operations
  reactor_recycle.jl        # PFR-with-recycle unit operations
  flowsheet.jl              # Stream/mix/split_stream + tear-stream solver
  reactor_staged_injection.jl # LDPE autoclave/tubular (staged initiator)
test/
  runtests.jl              # unit tests (known limits + EOS residual checks)
notebooks/
  ThermoExplorer.jl        # interactive Pluto UI
```

## Scope & honesty about the data

The Sanchez-Lacombe characteristic parameters (`Tstar`, `Pstar`,
`rhostar`) shipped in `SPECIES_DB` are typical literature-compilation
values, rounded for convenience, meant to make the package runnable out
of the box. They are **not** validated regression parameters for a
specific resin grade or solvent lot — refit them against real PVT data
before relying on this for engineering decisions. All of the theory
(Flory-Huggins activities/critical point/spinodal, the Sanchez-Lacombe
EOS itself) is implemented directly from its well-established closed-form
equations, and the test suite checks the code against known analytic
limits (e.g. the athermal, equal-size Flory-Huggins limit; the spinodal
touching the binodal at the critical point; the EOS residual vanishing at
the solved density).

The Sanchez-Lacombe binary mixing rules (`sl_mixing_rules`) were sourced
from a secondary review (Kontogeorgis, *A Survey of Equations of State
for Polymers*, IntechOpen 2012) and cross-checked by construction: they
reduce exactly to the pure-component parameters at the composition
limits, and the mixture's segment number, reconstructed from the mixed
parameters via the same formula used for pure components
(`segment_number`), reproduces the mixing rule's own value exactly — both
checked in the test suite.

Sanchez-Lacombe mixture chemical potentials/activities were deferred for
two prior increments of this package specifically because a literature
search turned up an explicit statement that published mixture
chemical-potential expressions for the Sanchez-Lacombe EOS are not always
thermodynamically consistent with each other. That statement traces to a
specific, citable result — von Konigslow, Park & Thompson, *Phys. Rev.
Applied* **8**, 044009 (2017) — which proves the inconsistency is
unavoidable for *any* composition-dependent hole-volume mixing rule
(including `sl_mixing_rules` above), but that chemical potentials *are*
thermodynamically consistent if the hole volume is held constant with
respect to composition when differentiating the free energy. This
package implements exactly that restricted, consistent case
(`sl_mixture_activities`, `src/sanchez_lacombe_activity.jl`), with its
own self-contained mixing rules (the source's Eqs. (18), (23), (36)) —
deliberately *not* reusing `sl_mixing_rules`, because an early attempt to
reuse it looked plausible (it passed a pure-component-limit check, which
almost any mixing rule satisfies) but failed a much stronger test: the
Euler relation `Σₖ nₖ μₖ − PV = F`, which must hold exactly whenever `μₖ`
and `P` are genuinely partial derivatives of the same free energy. The
chemical potential itself was independently rederived from the source's
free-energy expression by direct differentiation (matching the source's
own published formula), and the final implementation was verified against
that Euler relation numerically at several different states (not just
derived by hand) — both checked in the test suite. The trade-off, stated
plainly by the source and preserved here: this chemical potential does
not correspond to the *same* hole-volume convention as `sl_mixing_rules`'s
PVT predictions, so the two modules' outputs (density vs. activity) are
independently self-consistent but not mutually cross-consistent with a
single, unified equation of state.

## License

MIT — see [`LICENSE`](LICENSE).
