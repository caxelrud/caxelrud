"""
Sanchez-Lacombe binary mixture PVT behavior, via van der Waals-type mixing
rules for the characteristic parameters applied to the same reduced EOS
already implemented for pure components (`sanchez_lacombe.jl`).

Sources: the mixing rules below (segment/"hard-core" volume fractions `φᵢ`
from weight fractions, and quadratic combining rules for the characteristic
volume `ν*` and energy `ε*`, with binary interaction parameters `k12`,
`l12`) are as tabulated in Kontogeorgis's *A Survey of Equations of State
for Polymers* (IntechOpen, 2012), which cites them to Sanchez & Lacombe's
own mixture extension of their pure-fluid theory. The pure-component
reduced EOS this module builds on was independently cross-checked against
Vogiatzis, Megariotis & Theodorou, *J. Phys. Chem. B* (arXiv:1703.08983),
Appendix B, whose eq. B1 (Gibbs free energy) differentiates to exactly the
EOS already implemented in `sanchez_lacombe.jl` — a useful confirmation
that the underlying pure-fluid theory here is the standard one.

Scope note — chemical potentials/activities for Sanchez-Lacombe mixtures
are deliberately NOT implemented here. A literature search specifically
for this turned up an explicit warning that published mixture
chemical-potential expressions for the Sanchez-Lacombe EOS are not always
thermodynamically consistent with each other (Ferreira et al., *Fluid
Phase Equilib.* 2004 area of literature). Rather than pick one such
formula on faith, this module sticks to what can be verified directly:
mixture PVT behavior, which reduces to the correct pure-component limits
by construction (checked in the test suite) and reuses the
already-tested pure-component EOS solver rather than duplicating it.
"""

const _SL_MIX_R = 8.314462618  # J/(mol K)

_sl_segment_volume(sp::Species) = _SL_MIX_R * sp.Tstar / (sp.Pstar * 1e6)  # m^3/mol segments
_sl_segment_energy(sp::Species) = _SL_MIX_R * sp.Tstar                     # J/mol segments
_sl_segment_mass(sp::Species) = (sp.M / 1000) / segment_number(sp)         # kg/mol segments

"""
    sl_mixing_rules(sp1::Species, w1, sp2::Species, w2; k12=0.0, l12=0.0)

Sanchez-Lacombe binary mixing rules for the characteristic parameters.
`w1`, `w2` are weight (mass) fractions, `w1 + w2 == 1`. The segment
("hard-core volume") fractions are

``\\varphi_i = \\dfrac{w_i \\rho_i^* v_i^*}{\\sum_j w_j \\rho_j^* v_j^*}``

and the combining rules, with binary interaction parameters `k12`
(energy) and `l12` (volume):

``v_{12}^* = \\tfrac12(v_1^* + v_2^*)(1-l_{12}), \\qquad \\varepsilon_{12}^* = \\sqrt{\\varepsilon_1^*\\varepsilon_2^*}(1-k_{12})``

``v_{mix}^* = \\sum_i\\sum_j \\varphi_i\\varphi_j v_{ij}^*, \\qquad \\varepsilon_{mix}^* = \\dfrac{1}{v_{mix}^*}\\sum_i\\sum_j \\varphi_i\\varphi_j \\varepsilon_{ij}^* v_{ij}^*, \\qquad \\dfrac1{r_{mix}} = \\sum_i \\dfrac{\\varphi_i}{r_i}``

Returns a named tuple `(Tstar, Pstar, rhostar, M, r, φ1, φ2)` — the first
four using the same units/conventions as [`Species`](@ref) fields, `r` the
mixture segment number, and `φ1`/`φ2` the segment fractions. `M` is
constructed (`M = r_mix * Mseg_mix`) specifically so that
[`segment_number`](@ref) applied to a [`Species`](@ref) built from these
four fields reproduces `r` exactly — see [`sl_mixture_species`](@ref).
"""
function sl_mixing_rules(sp1::Species, w1, sp2::Species, w2; k12=0.0, l12=0.0)
    isapprox(w1 + w2, 1.0; atol=1e-8) || throw(ArgumentError("w1 + w2 must equal 1"))

    ν1, ν2 = _sl_segment_volume(sp1), _sl_segment_volume(sp2)
    ε1, ε2 = _sl_segment_energy(sp1), _sl_segment_energy(sp2)
    Mseg1, Mseg2 = _sl_segment_mass(sp1), _sl_segment_mass(sp2)
    r1, r2 = segment_number(sp1), segment_number(sp2)

    φ1 = w1 * Mseg1 / (w1 * Mseg1 + w2 * Mseg2)
    φ2 = 1 - φ1

    ν12 = 0.5 * (ν1 + ν2) * (1 - l12)
    ε12 = sqrt(ε1 * ε2) * (1 - k12)

    ν_mix = φ1^2 * ν1 + 2 * φ1 * φ2 * ν12 + φ2^2 * ν2
    ε_mix = (φ1^2 * ε1 * ν1 + 2 * φ1 * φ2 * ε12 * ν12 + φ2^2 * ε2 * ν2) / ν_mix
    r_mix = 1 / (φ1 / r1 + φ2 / r2)
    Mseg_mix = φ1 * Mseg1 + φ2 * Mseg2

    Tstar_mix = ε_mix / _SL_MIX_R
    Pstar_mix = ε_mix / ν_mix / 1e6   # Pa -> MPa
    rhostar_mix = Mseg_mix / ν_mix     # kg/m^3
    M_mix = r_mix * Mseg_mix * 1000    # kg/mol -> g/mol

    return (Tstar=Tstar_mix, Pstar=Pstar_mix, rhostar=rhostar_mix, M=M_mix, r=r_mix, φ1=φ1, φ2=φ2)
end

"""
    sl_mixture_species(sp1::Species, w1, sp2::Species, w2; k12=0.0, l12=0.0)

Build a "virtual" [`Species`](@ref) representing a binary Sanchez-Lacombe
mixture of `sp1` (weight fraction `w1`) and `sp2` (weight fraction
`w2 = 1 - w1`), via [`sl_mixing_rules`](@ref). The result can be fed
directly into any of the pure-component PVT functions in
`sanchez_lacombe.jl` (`reduced_density`, `density`, `specific_volume`,
`thermal_expansion_coefficient`, `isothermal_compressibility`) to get the
mixture's PVT behavior, reusing that already-tested EOS solver rather
than duplicating it.

`name`/`kind`/`δ`/`Vm` are placeholders (not meaningful for a mixture, and
not used by any PVT function); `Vm` is `NaN`.
"""
function sl_mixture_species(sp1::Species, w1, sp2::Species, w2; k12=0.0, l12=0.0)
    mix = sl_mixing_rules(sp1, w1, sp2, w2; k12=k12, l12=l12)
    return Species(
        name="$(sp1.name)/$(sp2.name) mixture (w1=$(round(w1, digits=3)))",
        kind=:solvent, M=mix.M, Vm=NaN, δ=NaN,
        Tstar=mix.Tstar, Pstar=mix.Pstar, rhostar=mix.rhostar,
    )
end
