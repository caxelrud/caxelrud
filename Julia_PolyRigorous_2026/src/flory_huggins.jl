"""
Flory-Huggins lattice theory for a binary polymer(2)-solvent(1) solution.

Convention: `φ₂` is the polymer volume fraction, `φ₁ = 1 - φ₂` the solvent
volume fraction, and `N` is the Flory-Huggins size ratio (number of
solvent-sized segments per polymer chain), i.e. `N = Vm(polymer)/Vm(solvent)`
— see [`degree_of_polymerization`](@ref).

All chemical-potential-type quantities below are the standard Flory (1953)
results:

``\\ln a_1 = \\ln\\phi_1 + \\left(1 - \\frac1N\\right)\\phi_2 + \\chi\\phi_2^2``

``\\ln a_2 = \\ln\\phi_2 - (N - 1)\\phi_1 + N\\chi\\phi_1^2``

where `a1` is the solvent activity and `a2` is the activity of the polymer
*per mole of chain* (not per segment).
"""

const R_GAS = 8.314462618  # J/(mol K)

"""
    chi_from_solubility(Vref, δ₁, δ₂, T; χₛ=0.34)

Estimate the Flory-Huggins interaction parameter from Hildebrand solubility
parameters (MPa^0.5), using the regular-solution enthalpic term plus an
empirical entropic correction `χₛ` (0.3-0.4 is a common rule-of-thumb range
for nonpolar polymer-solvent pairs).

`Vref` is the reference molar volume (commonly the solvent's), cm^3/mol.
`T` in K. Returns a dimensionless `χ`.
"""
function chi_from_solubility(Vref, δ₁, δ₂, T; χₛ=0.34)
    Vref_m3 = Vref * 1e-6      # cm^3/mol -> m^3/mol
    Δδ2 = (δ₁ - δ₂)^2 * 1e6    # (MPa^0.5)^2 -> Pa
    return Vref_m3 * Δδ2 / (R_GAS * T) + χₛ
end

"""
    ln_activity_solvent(φ₂, N, χ)

`ln(a1)` of the solvent, Flory-Huggins theory. `φ₂` is the polymer volume
fraction.
"""
function ln_activity_solvent(φ₂, N, χ)
    φ₁ = 1 - φ₂
    return log(φ₁) + (1 - 1 / N) * φ₂ + χ * φ₂^2
end

activity_solvent(φ₂, N, χ) = exp(ln_activity_solvent(φ₂, N, χ))

"""
    ln_activity_polymer(φ₂, N, χ)

`ln(a2)` of the polymer *per mole of chain*, Flory-Huggins theory.
"""
function ln_activity_polymer(φ₂, N, χ)
    φ₁ = 1 - φ₂
    return log(φ₂) - (N - 1) * φ₁ + N * χ * φ₁^2
end

activity_polymer(φ₂, N, χ) = exp(ln_activity_polymer(φ₂, N, χ))

"""
    flory_huggins_activities(φ₂, N, χ)

Convenience wrapper returning `(a1, a2)`.
"""
flory_huggins_activities(φ₂, N, χ) = (activity_solvent(φ₂, N, χ), activity_polymer(φ₂, N, χ))

"""
    osmotic_pressure(φ₂, N, χ, Vm1, T)

Solvent osmotic pressure `Π = -R T / Vm1 * ln(a1)` (Pa), with `Vm1` the
solvent molar volume in m^3/mol and `T` in K.
"""
function osmotic_pressure(φ₂, N, χ, Vm1, T)
    return -R_GAS * T / Vm1 * ln_activity_solvent(φ₂, N, χ)
end

"""
    chi_spinodal(φ₂, N)

Value of `χ` on the spinodal curve at polymer volume fraction `φ₂`, for
size ratio `N`, obtained from `d(ln a1)/d(φ₂) = 0`:

``\\chi_s(\\phi_2) = \\dfrac{\\frac{1}{1-\\phi_2} - \\left(1-\\frac1N\\right)}{2\\phi_2}``
"""
function chi_spinodal(φ₂, N)
    0 < φ₂ < 1 || throw(ArgumentError("φ₂ must be in (0, 1)"))
    return (1 / (1 - φ₂) - (1 - 1 / N)) / (2 * φ₂)
end

"""
    critical_point(N)

Flory-Huggins critical polymer volume fraction and interaction parameter for
a polymer(size `N`)-solvent(size 1) pair:

``\\phi_{2,c} = \\dfrac{1}{1+\\sqrt N}, \\qquad \\chi_c = \\dfrac12\\left(1+\\dfrac1{\\sqrt N}\\right)^2``

Returns a named tuple `(φ₂_c, χ_c)`.
"""
function critical_point(N)
    φ₂_c = 1 / (1 + sqrt(N))
    χ_c = 0.5 * (1 + 1 / sqrt(N))^2
    return (φ₂_c=φ₂_c, χ_c=χ_c)
end
