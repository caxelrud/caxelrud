"""
Liquid-liquid (polymer-solvent) phase-split calculations built on top of the
Flory-Huggins activity expressions in `flory_huggins.jl`.
"""

using NLsolve: nlsolve, converged

"""
    binodal_pair(N, chi; guess=nothing)

Solve for the two coexisting polymer volume fractions `(phi2a, phi2b)`,
`phi2a < phi2b`, at fixed Flory-Huggins `chi` and size ratio `N`, by imposing
equality of solvent and polymer activities between the two phases:

``a_1(\\phi_{2a}) = a_1(\\phi_{2b}), \\qquad a_2(\\phi_{2a}) = a_2(\\phi_{2b})``

Only meaningful for `chi` above the critical value (see [`critical_point`](@ref));
below it the only solution is the trivial `phi2a == phi2b`. When `guess` is
not given, an initial guess straddling the critical composition is
constructed automatically from `critical_point(N)`.

Returns a named tuple `(phi2a, phi2b, converged)`.
"""
function binodal_pair(N, chi; guess=nothing)
    if guess === nothing
        phi2c = critical_point(N).phi2c
        guess = (max(phi2c * 0.3, 1e-4), min(phi2c + (1 - phi2c) * 0.6, 1 - 1e-4))
    end
    function F!(Fvec, x)
        phi2a, phi2b = x
        if !(0 < phi2a < 1) || !(0 < phi2b < 1) || phi2a >= phi2b
            Fvec[1] = 1e6
            Fvec[2] = 1e6
            return Fvec
        end
        Fvec[1] = ln_activity_solvent(phi2a, N, chi) - ln_activity_solvent(phi2b, N, chi)
        Fvec[2] = ln_activity_polymer(phi2a, N, chi) - ln_activity_polymer(phi2b, N, chi)
        return Fvec
    end
    sol = nlsolve(F!, collect(Float64, guess); xtol=1e-13, ftol=1e-13, iterations=200)
    phi2a, phi2b = sol.zero
    if phi2a > phi2b
        phi2a, phi2b = phi2b, phi2a
    end
    return (phi2a=phi2a, phi2b=phi2b, converged=converged(sol))
end

"""
    binodal_curve(N; chi_max=nothing, npoints=40)

Trace the binodal (coexistence) curve for a polymer(size `N`)-solvent(size 1)
pair by continuation in `chi`, starting just above the critical point and
using each converged solution as the initial guess for the next. Returns a
`Vector` of named tuples `(chi, phi2a, phi2b)`, sorted by increasing `chi`.

`chi_max` defaults to `1.5 * chi_c`.
"""
function binodal_curve(N; chi_max=nothing, npoints=40)
    crit = critical_point(N)
    chic = crit.chic
    chi_hi = chi_max === nothing ? 1.5 * chic : chi_max
    chi_hi > chic || throw(ArgumentError("chi_max must exceed the critical chi ($chic)"))

    chis = range(chic * 1.001, chi_hi; length=npoints)
    guess = nothing

    out = NamedTuple{(:chi, :phi2a, :phi2b),Tuple{Float64,Float64,Float64}}[]
    for chi in chis
        res = binodal_pair(N, chi; guess=guess)
        res.converged || continue
        push!(out, (chi=chi, phi2a=res.phi2a, phi2b=res.phi2b))
        guess = (res.phi2a, res.phi2b)
    end
    return out
end

"""
    spinodal_curve(N; phi2_range=(1e-3, 0.999), npoints=200)

Sample the spinodal curve `chi_spinodal(phi2, N)` over `phi2_range`. Returns
a `Vector` of `(phi2, chi)` named tuples.
"""
function spinodal_curve(N; phi2_range=(1e-3, 0.999), npoints=200)
    phi2s = range(phi2_range[1], phi2_range[2]; length=npoints)
    return [(phi2=p, chi=chi_spinodal(p, N)) for p in phi2s]
end
