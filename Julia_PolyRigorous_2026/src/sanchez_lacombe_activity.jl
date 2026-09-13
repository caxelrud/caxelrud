"""
Sanchez-Lacombe binary mixture chemical potentials and activities.

`sanchez_lacombe_mixture.jl` deliberately stopped at mixture *PVT* behavior,
citing a literature warning that published Sanchez-Lacombe mixture
chemical-potential expressions are not always thermodynamically
consistent with each other. That warning traces to a specific, well
documented result: von Konigslow, Park & Thompson, *Phys. Rev. Applied*
**8**, 044009 (2017) show that the chemical potentials obtained from the
mixture free energy are thermodynamically inconsistent *for every*
composition-dependent hole-volume mixing rule — but *are* consistent for
a hole volume `v0` held constant with respect to composition when
differentiating the free energy (their Eq. (40) with its
mixing-rule-derivative term dropped). The trade-off, which they state
explicitly, is that a constant hole volume then does not exactly reduce
to each pure component's own hole volume at the composition limits —
*unless* the specific value chosen for `v0` already varies with
composition in a way that happens to hit the right pure-component value
at each limit, which is exactly what the mixing rule below (their own
Eqs. (18), (23), (36)) is built to do.

Important scope note: this module intentionally does **not** reuse
`sanchez_lacombe_mixture.jl`'s `sl_mixing_rules`/`sl_mixture_species`.
Those implement a *different*, independently cited (Kontogeorgis-review)
mixing-rule convention for `r`/`ε*`/`v0*` than the one this derivation
requires. An early version of this module tried to reuse them anyway and
looked plausible (it passed the pure-component-limit check, which nearly
any reasonable convention satisfies) — but failed a stronger check: the
Euler relation `Σκ nκ μκ - PV = F`, which must hold *exactly* for any
`μκ`, `P` that are genuinely `∂F/∂nκ` and `-∂F/∂V` of the *same* `F`. That
numerical check (reproduced in the test suite below) is what caught the
mismatch and is why this module computes its own self-contained mixing
rules from first principles instead.

Derivation and verification: rather than transcribe the source's Eq. (40),
we rederive the chemical potential directly, `μκ = (∂F/∂nκ)_{T,V,nκ'≠κ}`,
from the source's own free-energy expression (its Eq. (13)), holding `v0`
constant under the derivative (`∂v0/∂nκ = 0`). For component `κ` in a
binary mixture (`κ'` its partner):

``\\dfrac{\\mu_\\kappa}{RT} = -\\alpha_\\kappa(1 + \\ln\\phi_0) - 2\\alpha_\\kappa\\left(\\dfrac{\\epsilon_{\\kappa\\kappa}}{RT}\\phi_\\kappa + \\dfrac{\\epsilon_{\\kappa\\kappa'}}{RT}\\phi_{\\kappa'}\\right) + 1 + \\ln\\phi_\\kappa``

This independent rederivation reproduces the source's Eq. (40) exactly
once its mixing-rule-derivative term is dropped, and was additionally
checked (not just derived) by matching a finite-difference `∂F/∂nκ` of
the transcribed `F` against this closed form, and by the Euler-relation
check above — both reproduced as automated tests.

The quantities feeding that formula, all specific to this module (see the
docstrings below): segment-mole fractions `y_κ` and hard-core-volume
fractions `x_κ` built directly from mole numbers, segment counts
([`segment_number`](@ref)) and segment volumes (reusing the private
`_sl_segment_volume` helper already tested for pure components);
`v0 = y1 ν1 + y2 ν2` (Eq. (36)); `ε^*_{mix} = x_1^2\\epsilon_{11} +
2x_1x_2\\epsilon_{12} + x_2^2\\epsilon_{22}` (Eq. (18), with
`ε_{12} = \\sqrt{\\epsilon_{11}\\epsilon_{22}}(1-k_{12})`, the same
combining rule `sl_mixing_rules` uses); `T^*_{mix} = \\epsilon^*_{mix}/R`,
`P^*_{mix} = \\epsilon^*_{mix}/v_0/10^6`; and `1/r_{mix} = x_1/\\alpha_1 +
x_2/\\alpha_2` (Eq. (23), with `α_κ = r_κ ν_κ / v_0`). The mixture's
reduced density is then solved with the *same, already-tested* EOS
residual ([`sl_eos_residual`](@ref)) using these self-consistent
parameters, and `φ_κ = ρ̃ x_κ`, `φ_0 = 1 - ρ̃`.

`sl_mixture_activities` reports each component's activity *relative to
its own pure fluid at the same T, P* (solved with the already-tested
[`reduced_density`](@ref)), which is the practically useful quantity for
phase-equilibrium/solubility work — analogous to
[`flory_huggins_activities`](@ref) for the Flory-Huggins model. At the
pure-component composition limit (`w1 = 1` or `w1 = 0`) every mixing-rule
quantity above reduces exactly to its pure-component value by
construction, so `ln(a) = 0` there exactly — checked in the test suite.
"""

using Roots: find_zero, Bisection

# Common formula: μ_κ/(RT) = -α_κ(1+ln φ0) - 2α_κ(ε_κκ/RT φ_κ + ε_κκ'/RT φ_κ') + 1 + ln φ_κ
_sl_chem_pot_over_RT(α_κ, φ_κ, φ0, ε_self_over_RT, ε_cross_over_RT, φ_partner) =
    -α_κ * (1 + log(φ0)) - 2 * α_κ * (ε_self_over_RT * φ_κ + ε_cross_over_RT * φ_partner) + 1 + log(φ_κ)

# Pure-fluid reduction (v0 = the species' own segment volume, so α_κ = r_κ
# exactly, φ_κ = ρ̃, φ0 = 1-ρ̃, and there is no partner term).
function _sl_pure_chem_pot_over_RT(T, P, sp::Species)
    ρ̃ = reduced_density(T, P, sp)
    return _sl_chem_pot_over_RT(segment_number(sp), ρ̃, 1 - ρ̃, sp.Tstar / T, 0.0, 0.0)
end

"""
    sl_consistent_mixing_rules(sp1::Species, w1, sp2::Species, w2; k12=0.0)

The specific binary mixing rules needed for a thermodynamically consistent
Sanchez-Lacombe chemical potential (see the module docs): segment-mole
fractions `y1`/`y2`, hard-core-volume fractions `x1`/`x2`, the mixture's
segment/hole volume `v0` (Eq. (36), a segment-mole-fraction-weighted
average of the pure segment volumes), and the resulting `Tstar`, `Pstar`,
`r` (Eqs. (18), (23)). These are a *different* convention from
[`sl_mixing_rules`](@ref) — see the module docs for why they cannot be
mixed.

Returns a named tuple `(Tstar, Pstar, r, v0, x1, x2, y1, y2)`.
"""
function sl_consistent_mixing_rules(sp1::Species, w1, sp2::Species, w2; k12=0.0)
    isapprox(w1 + w2, 1.0; atol=1e-8) || throw(ArgumentError("w1 + w2 must equal 1"))

    r1, r2 = segment_number(sp1), segment_number(sp2)
    ν1, ν2 = _sl_segment_volume(sp1), _sl_segment_volume(sp2)

    n1_prop, n2_prop = w1 / sp1.M, w2 / sp2.M  # proportional to moles per unit mass
    y1 = (n1_prop * r1) / (n1_prop * r1 + n2_prop * r2)
    y2 = 1 - y1

    t1, t2 = n1_prop * r1 * ν1, n2_prop * r2 * ν2
    x1 = t1 / (t1 + t2)
    x2 = 1 - x1

    v0 = y1 * ν1 + y2 * ν2

    ε11 = _sl_segment_energy(sp1)
    ε22 = _sl_segment_energy(sp2)
    ε12 = sqrt(ε11 * ε22) * (1 - k12)
    εstar = x1^2 * ε11 + 2 * x1 * x2 * ε12 + x2^2 * ε22

    Tstar = εstar / _SL_MIX_R
    Pstar = εstar / v0 / 1e6

    α1, α2 = r1 * ν1 / v0, r2 * ν2 / v0
    r_mix = 1 / (x1 / α1 + x2 / α2)

    return (Tstar=Tstar, Pstar=Pstar, r=r_mix, v0=v0, x1=x1, x2=x2, y1=y1, y2=y2)
end

# Solve the SL EOS with an explicit r (not via a Species/segment_number),
# since sl_consistent_mixing_rules's r has no corresponding Species. Mirrors
# reduced_density's scan-then-bisect strategy exactly.
function _sl_reduced_density_explicit_r(T, P, Tstar, Pstar, r; ρ̃_bracket=(1e-6, 1 - 1e-9), nscan=400)
    T̃, P̃ = T / Tstar, P / Pstar
    f(ρ̃) = sl_eos_residual(ρ̃, T̃, P̃, r)
    lo, hi = ρ̃_bracket
    grid = range(hi, lo; length=nscan)
    f_prev = f(grid[1])
    for i in 2:length(grid)
        f_curr = f(grid[i])
        if f_prev * f_curr < 0
            return find_zero(f, (grid[i], grid[i-1]), Bisection())
        end
        f_prev = f_curr
    end
    throw(ErrorException("No liquid-branch root found for the mixture EOS residual at T=$T K, P=$P MPa."))
end

"""
    sl_mixture_ln_activities(T, P, sp1::Species, w1, sp2::Species, w2; k12=0.0)

`(ln_a1, ln_a2)`: the natural log of each component's activity in the
binary Sanchez-Lacombe mixture at temperature `T` (K) and pressure `P`
(MPa), relative to its own pure fluid at the same `T`, `P`. See the module
docs for the formula and its derivation/verification.
"""
function sl_mixture_ln_activities(T, P, sp1::Species, w1, sp2::Species, w2; k12=0.0)
    mix = sl_consistent_mixing_rules(sp1, w1, sp2, w2; k12=k12)

    ρ̃_mix = _sl_reduced_density_explicit_r(T, P, mix.Tstar, mix.Pstar, mix.r)
    φ1_actual = ρ̃_mix * mix.x1
    φ2_actual = ρ̃_mix * mix.x2
    φ0 = 1 - ρ̃_mix

    r1, r2 = segment_number(sp1), segment_number(sp2)
    ν1, ν2 = _sl_segment_volume(sp1), _sl_segment_volume(sp2)
    α1, α2 = r1 * ν1 / mix.v0, r2 * ν2 / mix.v0

    ε11_over_RT = sp1.Tstar / T
    ε22_over_RT = sp2.Tstar / T
    ε12_over_RT = sqrt(sp1.Tstar * sp2.Tstar) * (1 - k12) / T

    μ1_over_RT = _sl_chem_pot_over_RT(α1, φ1_actual, φ0, ε11_over_RT, ε12_over_RT, φ2_actual)
    μ2_over_RT = _sl_chem_pot_over_RT(α2, φ2_actual, φ0, ε22_over_RT, ε12_over_RT, φ1_actual)

    μ1_pure_over_RT = _sl_pure_chem_pot_over_RT(T, P, sp1)
    μ2_pure_over_RT = _sl_pure_chem_pot_over_RT(T, P, sp2)

    return (ln_a1=μ1_over_RT - μ1_pure_over_RT, ln_a2=μ2_over_RT - μ2_pure_over_RT)
end

"""
    sl_mixture_activities(T, P, sp1::Species, w1, sp2::Species, w2; k12=0.0)

`(a1, a2)`, i.e. `exp.(sl_mixture_ln_activities(...))`.
"""
function sl_mixture_activities(T, P, sp1::Species, w1, sp2::Species, w2; k12=0.0)
    r = sl_mixture_ln_activities(T, P, sp1, w1, sp2, w2; k12=k12)
    return (a1=exp(r.ln_a1), a2=exp(r.ln_a2))
end
