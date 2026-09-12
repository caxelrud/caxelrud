# PolyRigorous

A Julia + [Pluto.jl](https://plutojl.org/) toolkit for rigorous polymer
solution thermodynamics — the kind of calculations underlying commercial
tools like Aspen Polymers Plus. It's a from-scratch, open implementation,
not a clone of any proprietary source or model — see [Scope & honesty
about the data](#scope--honesty-about-the-data) below.

## What's here (v0.1.0)

This first version is a **thermodynamics core**:

- A small built-in species database (`species`, `SPECIES_DB`) with a
  handful of common polymers (polystyrene, PMMA, polyethylene) and
  solvents (toluene, benzene, cyclohexane, n-hexane).
- **Flory-Huggins lattice theory**: solvent/polymer activities, osmotic
  pressure, the spinodal curve, and the critical point, for a binary
  polymer-solvent system (`src/flory_huggins.jl`).
- **Liquid-liquid phase equilibrium**: binodal (coexistence) curve
  computation by solving the equal-activity conditions between two phases
  (`src/phase_equilibrium.jl`).
- **Sanchez-Lacombe lattice-fluid equation of state** for pure-component
  PVT behavior: density, specific volume, thermal expansion coefficient,
  isothermal compressibility (`src/sanchez_lacombe.jl`).
- An interactive **Pluto notebook**, `notebooks/ThermoExplorer.jl`, that
  puts sliders and dropdowns on top of all of the above.

## Roadmap (not yet implemented)

- Polymerization kinetics (free-radical, step-growth, coordination) for
  molecular weight distributions.
- Reactor unit operations (CSTR / PFR / batch) built on those kinetics.
- Sanchez-Lacombe *mixture* thermodynamics — binary mixing rules and the
  resulting chemical potentials/activities (the current version only
  covers pure-component PVT).
- A full flowsheet solver connecting multiple unit operations.

## Getting started

You need a working [Julia](https://julialang.org/) installation
(1.9+recommended).

```bash
git clone <this repo>
cd Julia_PolyRigorous_2026
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

This has been run end-to-end (Julia 1.13, all 39 tests passing) as part of building this package.

### Using the package directly

```julia
using Pkg
Pkg.activate(".")
using PolyRigorous

poly = species("polystyrene")
solv = species("toluene")
N = degree_of_polymerization(poly, solv)

chi = chi_from_solubility(solv.Vm, solv.delta, poly.delta, 298.15)
crit = critical_point(N)
@show chi, crit
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

## Project layout

```
Project.toml            # PolyRigorous package manifest
src/
  PolyRigorous.jl        # module entry point
  components.jl           # Species type + built-in database
  flory_huggins.jl         # Flory-Huggins activities, spinodal, critical point
  phase_equilibrium.jl     # binodal curve via NLsolve
  sanchez_lacombe.jl       # Sanchez-Lacombe pure-component EOS
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

## License

MIT — see [`LICENSE`](LICENSE).
