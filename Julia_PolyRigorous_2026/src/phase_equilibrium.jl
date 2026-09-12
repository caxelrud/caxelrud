"""
Liquid-liquid (polymer-solvent) phase-split calculations built on top of the
Flory-Huggins activity expressions in `flory_huggins.jl`.
"""

using NLsolve: nlsolve, converged
using Roots: find_zero, Bisection

"""
    _binodal_residual!(Fvec, x, N, chi)

Binodal residual directly in `(phi2a, phi2b)`.

Note this system is degenerate along the whole diagonal `phi2a == phi2b`
(both equal-activity conditions are then trivially `0 - 0 = 0`, for *any*
point on that line) — an unconstrained reparametrization that maps all of
`R^2` onto `(0,1)^2` (e.g. a logit/sigmoid transform) turns out to make
that trivial-solution manifold dominate the solver's basin of attraction,
pulling Newton onto it even from good starting points away from the
diagonal (this happened in practice — see the package tests). Working
directly in `(phi2a, phi2b)` and simply rejecting excursions outside
`(0, 1)` with a flat penalty avoids that failure mode; the guess (see
[`_spinodal_guess`](@ref)) is what has to do the work of staying in the
correct basin.
"""
function _binodal_residual!(Fvec, x, N, chi)
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

function _binodal_solve(N, chi, guess)
    sol = nlsolve(
        (Fvec, x) -> _binodal_residual!(Fvec, x, N, chi),
        collect(Float64, guess);
        xtol=1e-11, ftol=1e-11, iterations=400,
    )
    phi2a, phi2b = sol.zero
    if phi2a > phi2b
        phi2a, phi2b = phi2b, phi2a
    end
    return (phi2a=phi2a, phi2b=phi2b, converged=converged(sol))
end

"""
    _spinodal_guess(N, chi)

Construct a well-conditioned initial guess `(phi2a, phi2b)` for the binodal
at a given `chi > chi_c`, anchored to the spinodal curve rather than to a
fixed offset from the critical composition.

The spinodal `chi_spinodal(phi2, N)` (see `flory_huggins.jl`) is a *closed
form*, single-valued, monotonic function of `phi2` on each side of the
critical composition, so its two branches at a given `chi` are found by
plain 1D bisection — no delicate near-critical Newton continuation needed.
The true binodal always lies just outside the spinodal on each side, and
crucially the spinodal-branch separation itself shrinks to zero as
`chi -> chi_c` and grows away from it, which is exactly the scaling a
fixed-width guess was missing (see the package tests, which caught this:
a fixed-width guess converged to a spurious near-degenerate solution close
to the critical point).
"""
function _spinodal_guess(N, chi)
    phi2c = critical_point(N).phi2c
    f(phi2) = chi_spinodal(phi2, N) - chi
    phi2_lo = find_zero(f, (1e-10, phi2c - 1e-12), Bisection())
    phi2_hi = find_zero(f, (phi2c + 1e-12, 1 - 1e-10), Bisection())

    # Push outward from the spinodal by half the (signed) distance back to
    # the critical composition. This additive extrapolation is the accurate
    # one near chi_c, where phi2_lo/phi2_hi are themselves close to phi2c —
    # but it overshoots past 0 (or past 1 on the high side) once chi is well
    # above chi_c and phi2_lo/phi2_hi have shrunk toward the domain edges.
    # In that regime fall back to a fraction of phi2_lo itself (resp.
    # 1 - phi2_hi), which is exactly the scale-invariant behavior needed
    # there. Each formula is accurate exactly where the other breaks down.
    guess_lo_additive = 1.5 * phi2_lo - 0.5 * phi2c
    guess_lo = guess_lo_additive > 0.1 * phi2_lo ? guess_lo_additive : 0.5 * phi2_lo

    guess_hi_additive = 1.5 * phi2_hi - 0.5 * phi2c
    guess_hi = guess_hi_additive < 1 - 0.1 * (1 - phi2_hi) ? guess_hi_additive : 1 - 0.5 * (1 - phi2_hi)

    return (guess_lo, guess_hi)
end

"""
    binodal_pair(N, chi; guess=nothing)

Solve for the two coexisting polymer volume fractions `(phi2a, phi2b)`,
`phi2a < phi2b`, at fixed Flory-Huggins `chi` and size ratio `N`, by imposing
equality of solvent and polymer activities between the two phases:

``a_1(\\phi_{2a}) = a_1(\\phi_{2b}), \\qquad a_2(\\phi_{2a}) = a_2(\\phi_{2b})``

Only meaningful for `chi` above the critical value (see [`critical_point`](@ref));
below it the only solution is the trivial `phi2a == phi2b`. When `guess` is
not given, it is constructed automatically via [`_spinodal_guess`](@ref).

Returns a named tuple `(phi2a, phi2b, converged)`.
"""
function binodal_pair(N, chi; guess=nothing)
    crit = critical_point(N)
    chi > crit.chic || throw(ArgumentError("chi must exceed the critical chi ($(crit.chic)) for a nontrivial binodal to exist"))
    guess === nothing && (guess = _spinodal_guess(N, chi))
    return _binodal_solve(N, chi, guess)
end

"""
    binodal_curve(N; chi_max=nothing, npoints=40)

Trace the binodal (coexistence) curve for a polymer(size `N`)-solvent(size 1)
pair by solving [`binodal_pair`](@ref) independently at a grid of `chi`
values just above the critical point up to `chi_max` (each point uses its
own spinodal-anchored guess, so points don't depend on one another
converging). Returns a `Vector` of named tuples `(chi, phi2a, phi2b)`,
sorted by increasing `chi`.

`chi_max` defaults to `1.5 * chi_c`.
"""
function binodal_curve(N; chi_max=nothing, npoints=40)
    crit = critical_point(N)
    chic = crit.chic
    chi_hi = chi_max === nothing ? 1.5 * chic : chi_max
    chi_hi > chic || throw(ArgumentError("chi_max must exceed the critical chi ($chic)"))

    chis = range(chic * 1.001, chi_hi; length=npoints)

    out = NamedTuple{(:chi, :phi2a, :phi2b),Tuple{Float64,Float64,Float64}}[]
    for chi in chis
        res = binodal_pair(N, chi)
        res.converged || continue
        push!(out, (chi=chi, phi2a=res.phi2a, phi2b=res.phi2b))
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
