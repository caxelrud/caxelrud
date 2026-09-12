"""
Liquid-liquid (polymer-solvent) phase-split calculations built on top of the
Flory-Huggins activity expressions in `flory_huggins.jl`.
"""

using NLsolve: nlsolve, converged
using Roots: find_zero, Bisection

"""
    _binodal_residual!(Fvec, x, N, χ)

Binodal residual directly in `(φ₂a, φ₂b)`.

Note this system is degenerate along the whole diagonal `φ₂a == φ₂b` (both
equal-activity conditions are then trivially `0 - 0 = 0`, for *any* point
on that line) — an unconstrained reparametrization that maps all of `R^2`
onto `(0,1)^2` (e.g. a logit/sigmoid transform) turns out to make that
trivial-solution manifold dominate the solver's basin of attraction,
pulling Newton onto it even from good starting points away from the
diagonal (this happened in practice — see the package tests). Working
directly in `(φ₂a, φ₂b)` and simply rejecting excursions outside `(0, 1)`
with a flat penalty avoids that failure mode; the guess (see
[`_spinodal_guess`](@ref)) is what has to do the work of staying in the
correct basin.
"""
function _binodal_residual!(Fvec, x, N, χ)
    φ₂a, φ₂b = x
    if !(0 < φ₂a < 1) || !(0 < φ₂b < 1) || φ₂a >= φ₂b
        Fvec[1] = 1e6
        Fvec[2] = 1e6
        return Fvec
    end
    Fvec[1] = ln_activity_solvent(φ₂a, N, χ) - ln_activity_solvent(φ₂b, N, χ)
    Fvec[2] = ln_activity_polymer(φ₂a, N, χ) - ln_activity_polymer(φ₂b, N, χ)
    return Fvec
end

function _binodal_solve(N, χ, guess)
    sol = nlsolve(
        (Fvec, x) -> _binodal_residual!(Fvec, x, N, χ),
        collect(Float64, guess);
        xtol=1e-11, ftol=1e-11, iterations=400,
    )
    φ₂a, φ₂b = sol.zero
    if φ₂a > φ₂b
        φ₂a, φ₂b = φ₂b, φ₂a
    end
    return (φ₂a=φ₂a, φ₂b=φ₂b, converged=converged(sol))
end

"""
    _spinodal_guess(N, χ)

Construct a well-conditioned initial guess `(φ₂a, φ₂b)` for the binodal at
a given `χ > χ_c`, anchored to the spinodal curve rather than to a fixed
offset from the critical composition.

The spinodal `chi_spinodal(φ₂, N)` (see `flory_huggins.jl`) is a *closed
form*, single-valued, monotonic function of `φ₂` on each side of the
critical composition, so its two branches at a given `χ` are found by
plain 1D bisection — no delicate near-critical Newton continuation needed.
The true binodal always lies just outside the spinodal on each side, and
crucially the spinodal-branch separation itself shrinks to zero as
`χ -> χ_c` and grows away from it, which is exactly the scaling a
fixed-width guess was missing (see the package tests, which caught this:
a fixed-width guess converged to a spurious near-degenerate solution close
to the critical point).
"""
function _spinodal_guess(N, χ)
    φ₂_c = critical_point(N).φ₂_c
    f(φ₂) = chi_spinodal(φ₂, N) - χ
    φ₂_lo = find_zero(f, (1e-10, φ₂_c - 1e-12), Bisection())
    φ₂_hi = find_zero(f, (φ₂_c + 1e-12, 1 - 1e-10), Bisection())

    # Push outward from the spinodal by half the (signed) distance back to
    # the critical composition. This additive extrapolation is the accurate
    # one near χ_c, where φ₂_lo/φ₂_hi are themselves close to φ₂_c — but it
    # overshoots past 0 (or past 1 on the high side) once χ is well above
    # χ_c and φ₂_lo/φ₂_hi have shrunk toward the domain edges. In that
    # regime fall back to a fraction of φ₂_lo itself (resp. 1 - φ₂_hi),
    # which is exactly the scale-invariant behavior needed there. Each
    # formula is accurate exactly where the other breaks down.
    guess_lo_additive = 1.5 * φ₂_lo - 0.5 * φ₂_c
    guess_lo = guess_lo_additive > 0.1 * φ₂_lo ? guess_lo_additive : 0.5 * φ₂_lo

    guess_hi_additive = 1.5 * φ₂_hi - 0.5 * φ₂_c
    guess_hi = guess_hi_additive < 1 - 0.1 * (1 - φ₂_hi) ? guess_hi_additive : 1 - 0.5 * (1 - φ₂_hi)

    return (guess_lo, guess_hi)
end

"""
    binodal_pair(N, χ; guess=nothing)

Solve for the two coexisting polymer volume fractions `(φ₂a, φ₂b)`,
`φ₂a < φ₂b`, at fixed Flory-Huggins `χ` and size ratio `N`, by imposing
equality of solvent and polymer activities between the two phases:

``a_1(\\phi_{2a}) = a_1(\\phi_{2b}), \\qquad a_2(\\phi_{2a}) = a_2(\\phi_{2b})``

Only meaningful for `χ` above the critical value (see [`critical_point`](@ref));
below it the only solution is the trivial `φ₂a == φ₂b`. When `guess` is not
given, it is constructed automatically via [`_spinodal_guess`](@ref).

Returns a named tuple `(φ₂a, φ₂b, converged)`.
"""
function binodal_pair(N, χ; guess=nothing)
    crit = critical_point(N)
    χ > crit.χ_c || throw(ArgumentError("χ must exceed the critical χ ($(crit.χ_c)) for a nontrivial binodal to exist"))
    guess === nothing && (guess = _spinodal_guess(N, χ))
    return _binodal_solve(N, χ, guess)
end

"""
    binodal_curve(N; χ_max=nothing, npoints=40)

Trace the binodal (coexistence) curve for a polymer(size `N`)-solvent(size 1)
pair by solving [`binodal_pair`](@ref) independently at a grid of `χ`
values just above the critical point up to `χ_max` (each point uses its
own spinodal-anchored guess, so points don't depend on one another
converging). Returns a `Vector` of named tuples `(χ, φ₂a, φ₂b)`, sorted by
increasing `χ`.

`χ_max` defaults to `1.5 * χ_c`.
"""
function binodal_curve(N; χ_max=nothing, npoints=40)
    crit = critical_point(N)
    χ_c = crit.χ_c
    χ_hi = χ_max === nothing ? 1.5 * χ_c : χ_max
    χ_hi > χ_c || throw(ArgumentError("χ_max must exceed the critical χ ($χ_c)"))

    χs = range(χ_c * 1.001, χ_hi; length=npoints)

    out = NamedTuple{(:χ, :φ₂a, :φ₂b),Tuple{Float64,Float64,Float64}}[]
    for χ in χs
        res = binodal_pair(N, χ)
        res.converged || continue
        push!(out, (χ=χ, φ₂a=res.φ₂a, φ₂b=res.φ₂b))
    end
    return out
end

"""
    spinodal_curve(N; φ₂_range=(1e-3, 0.999), npoints=200)

Sample the spinodal curve `chi_spinodal(φ₂, N)` over `φ₂_range`. Returns a
`Vector` of `(φ₂, χ)` named tuples.
"""
function spinodal_curve(N; φ₂_range=(1e-3, 0.999), npoints=200)
    φ₂s = range(φ₂_range[1], φ₂_range[2]; length=npoints)
    return [(φ₂=p, χ=chi_spinodal(p, N)) for p in φ₂s]
end
