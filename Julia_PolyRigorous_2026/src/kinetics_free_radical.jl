"""
Free-radical polymerization kinetics: initiation, propagation, termination,
and the resulting molecular weight averages, under the standard
quasi-steady-state approximation (QSSA) for the radical concentration.

Convention follows Odian, *Principles of Polymerization*: initiator `I`,
monomer `M`, initiator decomposition rate constant `kd`, initiator
efficiency `f` (fraction of primary radicals that successfully start a
chain), propagation rate constant `kp`, termination rate constant
`kt = ktc + ktd` (combination + disproportionation).
"""

"""
    initiation_rate(f, kd, I)

Rate of radical generation `Ri = 2 f kd [I]` (mol/(L·s), for concentrations
in mol/L and `kd` in 1/s).
"""
initiation_rate(f, kd, I) = 2 * f * kd * I

"""
    radical_concentration(f, kd, I, kt)

Steady-state radical concentration `[M.] = sqrt(f kd [I] / kt)`, obtained
from the QSSA `Ri = Rt = 2 kt [M.]^2`.
"""
radical_concentration(f, kd, I, kt) = sqrt(f * kd * I / kt)

"""
    propagation_rate(kp, M, Mrad)

Rate of monomer consumption by propagation, `Rp = kp [M] [M.]`.
"""
propagation_rate(kp, M, Mrad) = kp * M * Mrad

"""
    kinetic_chain_length(kp, M, f, kd, I, kt)

Kinetic chain length `nu = Rp / Ri = kp [M] / (2 sqrt(kt f kd [I]))`: the
average number of monomer units consumed per radical that initiates a
chain.
"""
kinetic_chain_length(kp, M, f, kd, I, kt) = kp * M / (2 * sqrt(kt * f * kd * I))

"""
    Xn_combination(nu)

Number-average degree of polymerization when all termination is by
combination: `Xn = 2 nu` (two kinetic chains join into one dead chain).
"""
Xn_combination(nu) = 2 * nu

"""
    Xn_disproportionation(nu)

Number-average degree of polymerization when all termination is by
disproportionation: `Xn = nu` (each kinetic chain becomes its own dead
chain).
"""
Xn_disproportionation(nu) = nu

"""
    Xn_mixed(nu, delta)

Number-average degree of polymerization for termination split between
disproportionation (fraction `delta`, `0 <= delta <= 1`) and combination
(fraction `1 - delta`): `Xn = 2 nu / (1 + delta)`. Reduces to
[`Xn_combination`](@ref) at `delta = 0` and [`Xn_disproportionation`](@ref)
at `delta = 1`.
"""
function Xn_mixed(nu, delta)
    0 <= delta <= 1 || throw(ArgumentError("delta must be in [0, 1]"))
    return 2 * nu / (1 + delta)
end

"""
    monomer_concentration(t, kp, kd, kt, f, I0, M0)

Monomer concentration at time `t` in an isothermal batch reactor, from the
closed-form solution of `d[M]/dt = -kp [M] [M.]` under QSSA with a slowly
decaying initiator `[I](t) = I0 exp(-kd t)`:

``\\ln\\frac{[M]_0}{[M](t)} = \\frac{2 k_p}{k_d}\\sqrt{f k_d [I]_0 / k_t}\\,\\bigl(1 - e^{-k_d t/2}\\bigr)``

(Odian, *Principles of Polymerization*).
"""
function monomer_concentration(t, kp, kd, kt, f, I0, M0)
    ln_ratio = (2 * kp / kd) * sqrt(f * kd * I0 / kt) * (1 - exp(-kd * t / 2))
    return M0 * exp(-ln_ratio)
end

"""
    conversion(t, kp, kd, kt, f, I0, M0)

Fractional monomer conversion `1 - [M](t)/[M]_0` at time `t`, from
[`monomer_concentration`](@ref).
"""
conversion(t, kp, kd, kt, f, I0, M0) = 1 - monomer_concentration(t, kp, kd, kt, f, I0, M0) / M0
