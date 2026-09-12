"""
Sanchez-Lacombe lattice-fluid equation of state (pure component).

Sanchez & Lacombe, J. Phys. Chem. 80, 2352 (1976). Reduced variables are
defined per-species from the characteristic parameters `Tstar`, `Pstar`,
`rhostar` (see [`Species`](@ref)):

``\\tilde T = T/T^*, \\qquad \\tilde P = P/P^*, \\qquad \\tilde\\rho = \\rho/\\rho^*``

and the equation of state is

``\\tilde\\rho^2 + \\tilde P + \\tilde T\\left[\\ln(1-\\tilde\\rho) + \\left(1-\\dfrac1r\\right)\\tilde\\rho\\right] = 0``

where `r` is the number of lattice sites per molecule ([`segment_number`](@ref)).
The code below spells these reduced variables exactly as in the math,
`T̃`, `P̃`, `ρ̃` (a base letter plus a combining tilde, which Julia accepts
directly in identifiers).

Scope note: this module covers *pure-component* PVT behavior only (reduced
density, specific volume, thermal expansion, isothermal compressibility),
all obtained directly from the EOS above. Mixture (polymer + solvent)
Sanchez-Lacombe thermodynamics — mixing rules for `Tstar`/`Pstar`/`rhostar`
and the resulting chemical potentials/activities — is deliberately left for
a later version rather than shipped from an unverified closed-form
expression; see the README roadmap.
"""

using Roots: find_zero, Bisection

"""
    sl_eos_residual(ρ̃, T̃, P̃, r)

Left-hand side of the Sanchez-Lacombe EOS; zero at the physical reduced
density.
"""
function sl_eos_residual(ρ̃, T̃, P̃, r)
    return ρ̃^2 + P̃ + T̃ * (log1p(-ρ̃) + (1 - 1 / r) * ρ̃)
end

"""
    reduced_density(T, P, sp::Species; ρ̃_bracket=(1e-6, 1 - 1e-9), nscan=400)

Solve the Sanchez-Lacombe EOS for the liquid/melt-branch reduced density
`ρ̃ = ρ/rhostar` of species `sp` at temperature `T` (K) and pressure `P`
(MPa).

Like other cubic-type equations of state (van der Waals, Peng-Robinson,
...), the Sanchez-Lacombe EOS residual is not monotonic in `ρ̃` over the
whole `(0, 1)` range in general: besides the dense liquid/melt root (close
to the close-packed limit `rhostar`), there can also be a low-density,
vapor-like root. This package targets condensed polymer melts and
polymer-solvent solutions, so we want the liquid branch specifically. To
get it reliably (rather than risk a naive bisection over the full range
latching onto the wrong root), we scan `ρ̃` down from the high-density end
and bisect within the *first* sign change encountered — i.e. the root
nearest the close-packed limit.
"""
function reduced_density(T, P, sp::Species; ρ̃_bracket=(1e-6, 1 - 1e-9), nscan=400)
    r = segment_number(sp)
    T̃ = T / sp.Tstar
    P̃ = P / sp.Pstar
    f(ρ̃) = sl_eos_residual(ρ̃, T̃, P̃, r)

    lo, hi = ρ̃_bracket
    grid = range(hi, lo; length=nscan)  # scan from close-packed limit downward
    f_prev = f(grid[1])
    for i in 2:length(grid)
        f_curr = f(grid[i])
        if f_prev * f_curr < 0
            return find_zero(f, (grid[i], grid[i-1]), Bisection())
        end
        f_prev = f_curr
    end
    throw(ErrorException(
        "No liquid-branch root found for the EOS residual on ρ̃ in ($lo, $hi) at " *
        "T=$T K, P=$P MPa for $(sp.name); the state point may be outside this solver's " *
        "validity range."
    ))
end

"""
    density(T, P, sp::Species)

Mass density (kg/m^3) of species `sp` at `T` (K), `P` (MPa) from the
Sanchez-Lacombe EOS.
"""
density(T, P, sp::Species) = reduced_density(T, P, sp) * sp.rhostar

"""
    specific_volume(T, P, sp::Species)

Specific volume (m^3/kg), i.e. `1/density`.
"""
specific_volume(T, P, sp::Species) = 1 / density(T, P, sp)

"""
    thermal_expansion_coefficient(T, P, sp::Species; dT=1e-2)

Isobaric thermal expansion coefficient `α = (1/V) dV/dT` (1/K), by central
finite difference on [`specific_volume`](@ref).
"""
function thermal_expansion_coefficient(T, P, sp::Species; dT=1e-2)
    Vp = specific_volume(T + dT, P, sp)
    Vm = specific_volume(T - dT, P, sp)
    V0 = specific_volume(T, P, sp)
    return (Vp - Vm) / (2dT) / V0
end

"""
    isothermal_compressibility(T, P, sp::Species; dP=1e-3)

Isothermal compressibility `κ = -(1/V) dV/dP` (1/MPa), by central finite
difference on [`specific_volume`](@ref). `dP` in MPa; requires `P > dP`.
"""
function isothermal_compressibility(T, P, sp::Species; dP=1e-3)
    P > dP || throw(ArgumentError("P must exceed dP for a centered finite difference"))
    Vp = specific_volume(T, P + dP, sp)
    Vm = specific_volume(T, P - dP, sp)
    V0 = specific_volume(T, P, sp)
    return -(Vp - Vm) / (2dP) / V0
end
