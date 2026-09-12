"""
Reactor unit operations built on top of `kinetics_free_radical.jl` and
`kinetics_step_growth.jl`: ideal, isothermal, steady single CSTR and PFR.

Convention: `τ` is the reactor residence time (space time) `V/Q` — reactor
volume over volumetric feed flow rate. Subscript `_in` marks a feed
(inlet) concentration, as opposed to the reactor/outlet value.

An ideal PFR is kinetically *equivalent* to a batch reactor run for a time
equal to its residence time: with no back-mixing, each fluid element
travels down the PFR exactly like a closed batch reactor evolving from
`t=0` to `t=τ`. So the PFR functions below are thin wrappers around the
already-implemented batch kinetics, exposed under reactor-engineering
vocabulary (`τ`) for symmetry with the CSTR functions — not new physics.

A CSTR, by contrast, is a genuinely different (algebraic, not
differential) problem: the whole vessel sits at one steady-state
composition, fed continuously.
"""

using Roots: find_zero, Bisection

# ----------------------------------------------------------------------------
# Free-radical polymerization
# ----------------------------------------------------------------------------

"""
    cstr_free_radical(τ, k_p, k_d, k_t, f, I_in, M_in)

Steady-state composition of an ideal, isothermal CSTR running free-radical
polymerization, from the steady-state species balances
`(Cin - C)/τ = (consumption rate of C)`:

- Initiator: `(I_in - I)/τ = k_d I`, giving `I = I_in / (1 + k_d τ)`.
- Monomer: `(M_in - M)/τ = k_p M [M.]`, with `[M.]` the QSSA radical
  concentration ([`radical_concentration`](@ref)) evaluated at the
  reactor's own (outlet) initiator concentration `I` — valid because
  radical lifetimes (µs-ms) are always far shorter than reactor residence
  times, so radicals equilibrate to the *local* steady state. Since `[M.]`
  doesn't depend on `M`, this balance is linear in `M` and solves in
  closed form: `M = M_in / (1 + τ k_p [M.])`.

Returns a named tuple `(I, M, conversion)`.
"""
function cstr_free_radical(τ, k_p, k_d, k_t, f, I_in, M_in)
    I = I_in / (1 + k_d * τ)
    Mrad = radical_concentration(f, k_d, I, k_t)
    M = M_in / (1 + τ * k_p * Mrad)
    return (I=I, M=M, conversion=1 - M / M_in)
end

"""
    pfr_free_radical(τ, k_p, k_d, k_t, f, I_in, M_in)

Outlet composition of an ideal PFR with residence time `τ`, running
free-radical polymerization — identical to a batch reactor at reaction
time `t = τ` (see module docs); a thin wrapper around
[`monomer_concentration`](@ref).

Returns a named tuple `(M, conversion)`.
"""
function pfr_free_radical(τ, k_p, k_d, k_t, f, I_in, M_in)
    M = monomer_concentration(τ, k_p, k_d, k_t, f, I_in, M_in)
    return (M=M, conversion=1 - M / M_in)
end

# ----------------------------------------------------------------------------
# Step-growth polymerization
# ----------------------------------------------------------------------------

"""
    cstr_step_growth_external_catalyst(τ, k, c_in)

Steady-state functional-group concentration of an ideal, isothermal CSTR
running externally-catalyzed step-growth polymerization (second order
overall, see [`extent_reaction_external_catalyst`](@ref)), from
`(c_in - c)/τ = k c^2`, a quadratic in `c` with the physical (non-negative)
root

``c = \\dfrac{-1 + \\sqrt{1 + 4 k \\tau c_{in}}}{2 k \\tau} = \\dfrac{2 c_{in}}{1 + \\sqrt{1 + 4 k \\tau c_{in}}}``

using the rationalized (numerically stable as `τ -> 0`) form on the right
— the first form subtracts two nearly-equal numbers there and produces
`0/0`.

Returns a named tuple `(c, p)` where `p = 1 - c/c_in` is the extent of
reaction.
"""
function cstr_step_growth_external_catalyst(τ, k, c_in)
    c = 2 * c_in / (1 + sqrt(1 + 4 * k * τ * c_in))
    return (c=c, p=1 - c / c_in)
end

"""
    cstr_step_growth_self_catalyzed(τ, k, c_in)

Steady-state functional-group concentration of an ideal, isothermal CSTR
running self-catalyzed step-growth polymerization (third order overall,
see [`extent_reaction_self_catalyzed`](@ref)), from
`(c_in - c)/τ = k c^3`. Unlike the external-catalyst case this cubic has
no convenient closed form here, but `g(c) = k τ c^3 + c - c_in` is
strictly increasing in `c` (`g' = 3 k τ c^2 + 1 > 0` always), with
`g(0) = -c_in < 0` and `g(c_in) = k τ c_in^3 > 0`, so it has exactly one
root in `(0, c_in)`, found reliably by bisection.

Returns a named tuple `(c, p)` where `p = 1 - c/c_in` is the extent of
reaction.
"""
function cstr_step_growth_self_catalyzed(τ, k, c_in)
    g(c) = k * τ * c^3 + c - c_in
    c = find_zero(g, (0.0, c_in), Bisection())
    return (c=c, p=1 - c / c_in)
end

"""
    pfr_step_growth_external_catalyst(τ, k, c_in)

Outlet extent of reaction of an ideal PFR with residence time `τ`, running
externally-catalyzed step-growth polymerization — identical to a batch
reactor at reaction time `t = τ`; a thin wrapper around
[`extent_reaction_external_catalyst`](@ref).
"""
pfr_step_growth_external_catalyst(τ, k, c_in) = (p=extent_reaction_external_catalyst(k, c_in, τ),)

"""
    pfr_step_growth_self_catalyzed(τ, k, c_in)

Outlet extent of reaction of an ideal PFR with residence time `τ`, running
self-catalyzed step-growth polymerization — identical to a batch reactor
at reaction time `t = τ`; a thin wrapper around
[`extent_reaction_self_catalyzed`](@ref).
"""
pfr_step_growth_self_catalyzed(τ, k, c_in) = (p=extent_reaction_self_catalyzed(k, c_in, τ),)
