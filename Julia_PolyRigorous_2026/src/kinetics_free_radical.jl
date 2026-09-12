"""
Free-radical polymerization kinetics: initiation, propagation, termination,
and the resulting molecular weight averages, under the standard
quasi-steady-state approximation (QSSA) for the radical concentration.

Convention follows Odian, *Principles of Polymerization*: initiator `I`,
monomer `M`, initiator decomposition rate constant `k_d`, initiator
efficiency `f` (fraction of primary radicals that successfully start a
chain), propagation rate constant `k_p`, termination rate constant
`k_t = k_tc + k_td` (combination + disproportionation).
"""

"""
    initiation_rate(f, k_d, I)

Rate of radical generation `Ri = 2 f k_d [I]` (mol/(L·s), for concentrations
in mol/L and `k_d` in 1/s).
"""
initiation_rate(f, k_d, I) = 2 * f * k_d * I

"""
    radical_concentration(f, k_d, I, k_t)

Steady-state radical concentration `[M.] = sqrt(f k_d [I] / k_t)`, obtained
from the QSSA `Ri = Rt = 2 k_t [M.]^2`.
"""
radical_concentration(f, k_d, I, k_t) = sqrt(f * k_d * I / k_t)

"""
    propagation_rate(k_p, M, Mrad)

Rate of monomer consumption by propagation, `Rp = k_p [M] [M.]`.
"""
propagation_rate(k_p, M, Mrad) = k_p * M * Mrad

"""
    kinetic_chain_length(k_p, M, f, k_d, I, k_t)

Kinetic chain length `ν = Rp / Ri = k_p [M] / (2 sqrt(k_t f k_d [I]))`: the
average number of monomer units consumed per radical that initiates a
chain.
"""
kinetic_chain_length(k_p, M, f, k_d, I, k_t) = k_p * M / (2 * sqrt(k_t * f * k_d * I))

"""
    Xn_combination(ν)

Number-average degree of polymerization when all termination is by
combination: `Xn = 2 ν` (two kinetic chains join into one dead chain).
"""
Xn_combination(ν) = 2 * ν

"""
    Xn_disproportionation(ν)

Number-average degree of polymerization when all termination is by
disproportionation: `Xn = ν` (each kinetic chain becomes its own dead
chain).
"""
Xn_disproportionation(ν) = ν

"""
    Xn_mixed(ν, δ)

Number-average degree of polymerization for termination split between
disproportionation (fraction `δ`, `0 <= δ <= 1`) and combination (fraction
`1 - δ`): `Xn = 2 ν / (1 + δ)`. Reduces to [`Xn_combination`](@ref) at
`δ = 0` and [`Xn_disproportionation`](@ref) at `δ = 1`.
"""
function Xn_mixed(ν, δ)
    0 <= δ <= 1 || throw(ArgumentError("δ must be in [0, 1]"))
    return 2 * ν / (1 + δ)
end

"""
    monomer_concentration(t, k_p, k_d, k_t, f, I0, M0)

Monomer concentration at time `t` in an isothermal batch reactor, from the
closed-form solution of `d[M]/dt = -k_p [M] [M.]` under QSSA with a slowly
decaying initiator `[I](t) = I0 exp(-k_d t)`:

``\\ln\\frac{[M]_0}{[M](t)} = \\frac{2 k_p}{k_d}\\sqrt{f k_d [I]_0 / k_t}\\,\\bigl(1 - e^{-k_d t/2}\\bigr)``

(Odian, *Principles of Polymerization*).
"""
function monomer_concentration(t, k_p, k_d, k_t, f, I0, M0)
    ln_ratio = (2 * k_p / k_d) * sqrt(f * k_d * I0 / k_t) * (1 - exp(-k_d * t / 2))
    return M0 * exp(-ln_ratio)
end

"""
    conversion(t, k_p, k_d, k_t, f, I0, M0)

Fractional monomer conversion `1 - [M](t)/[M]_0` at time `t`, from
[`monomer_concentration`](@ref).
"""
conversion(t, k_p, k_d, k_t, f, I0, M0) = 1 - monomer_concentration(t, k_p, k_d, k_t, f, I0, M0) / M0
