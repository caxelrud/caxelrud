"""
PFR-with-recycle unit operations: a single ideal PFR whose outlet stream is
split, with a fraction recycled back and mixed with the fresh feed before
re-entering the reactor. This is the first step beyond `reactors.jl`'s
single-pass and reactor-train models toward genuine flowsheet topology: the
recycle stream's composition depends on the very reactor outlet it feeds,
so solving it is an implicit ("flowsheet-style") problem rather than a
single forward calculation.

Convention (Levenspiel, *Chemical Reaction Engineering*): `τ = V/Q_fresh`
is the space time based on the *fresh* feed flow rate alone (not the total
flow through the reactor), and `R = Q_recycle/Q_fresh` is the recycle
ratio. The reactor itself sees a *larger* total flow, `Q_fresh(1+R)`, so
its actual (internal) holding time is shorter:

``\\tau_{hold} = \\dfrac{\\tau}{1+R}``

A stream splitter doesn't change composition, so the recycle stream has
exactly the reactor's own outlet composition. A stream mixer combines the
fresh feed and the recycle stream by flow-weighted averaging:

``C_{in} = \\dfrac{C_{fresh} + R\\,C_{out}}{1+R}``

The reactor is an ideal PFR, so `C_out` is the batch/PFR-equivalent
evolution of `C_in` over `τ_hold` — but `C_in` itself depends on `C_out`
through the mixing equation above, so `C_out` appears on both sides. For
the free-radical and coordination cases below, the underlying batch
kinetics happen to be *linear* in the initial concentration, which turns
the recycle mass balance into an explicitly solvable equation for `C_out`.
For step-growth kinetics the batch kinetics are nonlinear, so the
external-catalyst case (still solvable in closed form, via a quadratic)
and the self-catalyzed case (solved numerically — no closed form found)
are handled differently — see each function's docstring.

Every function here reduces to the corresponding plain PFR (`pfr_*` in
`reactors.jl`) at `R = 0` by construction (`C_in = C_fresh`,
`τ_hold = τ`), and approaches the corresponding CSTR (`cstr_*`, same `τ`)
in the large-`R` limit — the classic reactor-engineering result that a PFR
with infinite recycle behaves like a CSTR of the same nominal residence
time — checked directly in the test suite rather than assumed.
"""

using Roots: find_zero, Bisection

# ----------------------------------------------------------------------------
# Free-radical polymerization
# ----------------------------------------------------------------------------

"""
    pfr_free_radical_recycle(τ, R, k_p, k_d, k_t, f, I_fresh, M_fresh)

PFR with recycle ratio `R` running free-radical polymerization. Both the
initiator and monomer balances are linear in the reactor's inlet
concentration (batch decay for the initiator, and
[`monomer_concentration`](@ref) is linear in its `M0` argument — see
module docs), so both solve in closed form.

Initiator: writing `α = exp(-k_d τ_hold)` for the batch decay factor over
the internal holding time `τ_hold = τ/(1+R)`, combining `I_out = α I_in`
with the mixing equation `I_in = (I_fresh + R I_out)/(1+R)` and solving for
`I_out` gives

``I_{out} = \\dfrac{I_{fresh}\\,\\alpha}{(1+R) - R\\alpha}``

Monomer: writing `κ` for the batch conversion factor
`monomer_concentration(τ_hold, k_p, k_d, k_t, f, I_in, 1.0)` (evaluated at
the reactor's own inlet initiator concentration `I_in`, found above),
`M_out = κ M_in` is linear in `M_in`, so combining with the mixing equation
the same way gives

``M_{out} = \\dfrac{M_{fresh}\\,\\kappa}{(1+R) - R\\kappa}``

Both denominators are `>= 1 > 0` (since `0 <= α, κ <= 1` and `R >= 0`), so
neither expression is ever singular, including at `R = 0` (where they
reduce to the plain batch-decay and [`pfr_free_radical`](@ref) formulas
exactly, since `α, κ` are then evaluated at `τ_hold = τ`).

Returns a named tuple `(I_in, I_out, M_in, M_out, conversion)`, with
`conversion = 1 - M_out/M_fresh` against the overall fresh feed.
"""
function pfr_free_radical_recycle(τ, R, k_p, k_d, k_t, f, I_fresh, M_fresh)
    τ_hold = τ / (1 + R)
    α = exp(-k_d * τ_hold)
    I_out = I_fresh * α / ((1 + R) - R * α)
    I_in = (I_fresh + R * I_out) / (1 + R)

    κ = monomer_concentration(τ_hold, k_p, k_d, k_t, f, I_in, 1.0)
    M_out = M_fresh * κ / ((1 + R) - R * κ)
    M_in = (M_fresh + R * M_out) / (1 + R)

    return (I_in=I_in, I_out=I_out, M_in=M_in, M_out=M_out, conversion=1 - M_out / M_fresh)
end

# ----------------------------------------------------------------------------
# Coordination polymerization
# ----------------------------------------------------------------------------

"""
    pfr_coordination_recycle(τ, R, k_p, k_t, C0_star_fresh, M_fresh)

PFR with recycle ratio `R` running coordination polymerization —
structurally identical to [`pfr_free_radical_recycle`](@ref), since active
sites decay by simple first-order kinetics (`k_t` here plays the role
`k_d` plays for the free-radical initiator) and
[`monomer_concentration_coordination`](@ref) is likewise linear in its
`M0` argument.

Returns a named tuple `(C_star_in, C_star_out, M_in, M_out, conversion)`,
with `conversion = 1 - M_out/M_fresh` against the overall fresh feed.
"""
function pfr_coordination_recycle(τ, R, k_p, k_t, C0_star_fresh, M_fresh)
    τ_hold = τ / (1 + R)
    α = exp(-k_t * τ_hold)
    C_star_out = C0_star_fresh * α / ((1 + R) - R * α)
    C_star_in = (C0_star_fresh + R * C_star_out) / (1 + R)

    κ = monomer_concentration_coordination(τ_hold, k_p, C_star_in, k_t, 1.0)
    M_out = M_fresh * κ / ((1 + R) - R * κ)
    M_in = (M_fresh + R * M_out) / (1 + R)

    return (C_star_in=C_star_in, C_star_out=C_star_out, M_in=M_in, M_out=M_out, conversion=1 - M_out / M_fresh)
end

# ----------------------------------------------------------------------------
# Step-growth polymerization
# ----------------------------------------------------------------------------

"""
    pfr_step_growth_external_catalyst_recycle(τ, R, k, c_fresh)

PFR with recycle ratio `R` running externally-catalyzed (second-order)
step-growth polymerization. Writing `κ = k τ_hold` and
`c_{in}(x) = (c_{fresh} + R x)/(1+R)` for the mixing equation, the recycle
condition `x = c_{in}(x)/(1 + κ c_{in}(x))` (the batch/PFR second-order
decay, [`extent_reaction_external_catalyst`](@ref), applied to the
reactor's own inlet) rearranges to a quadratic in `x`:

``R\\kappa\\,x^2 + (1 + c_{fresh}\\kappa)\\,x - c_{fresh} = 0``

with physical (non-negative) root, in the rationalized (numerically stable
as `R -> 0` or `κ -> 0`) form:

``x = \\dfrac{2 c_{fresh}}{(1+c_{fresh}\\kappa) + \\sqrt{(1+c_{fresh}\\kappa)^2 + 4 R \\kappa c_{fresh}}}``

which reduces to the plain PFR closed form at `R = 0` (the `Rκx²` term
vanishes, leaving the linear equation directly) without any
special-casing.

Returns a named tuple `(c_in, c_out, p)`, with `p = 1 - c_out/c_fresh` the
overall extent of reaction against the fresh feed.
"""
function pfr_step_growth_external_catalyst_recycle(τ, R, k, c_fresh)
    τ_hold = τ / (1 + R)
    κ = k * τ_hold
    c_out = 2 * c_fresh / ((1 + c_fresh * κ) + sqrt((1 + c_fresh * κ)^2 + 4 * R * κ * c_fresh))
    c_in = (c_fresh + R * c_out) / (1 + R)
    return (c_in=c_in, c_out=c_out, p=1 - c_out / c_fresh)
end

"""
    pfr_step_growth_self_catalyzed_recycle(τ, R, k, c_fresh)

PFR with recycle ratio `R` running self-catalyzed (third-order) step-growth
polymerization. Unlike the external-catalyst case, mixing the third-order
batch decay ([`extent_reaction_self_catalyzed`](@ref)) with the recycle
condition does not reduce to a low-order polynomial with a convenient
closed form, so `x = c_out` is found by root-finding on

``h(x) = x - \\dfrac{c_{in}(x)}{\\sqrt{1 + 2\\kappa\\,c_{in}(x)^2}}, \\qquad c_{in}(x) = \\dfrac{c_{fresh}+Rx}{1+R}``

`h` is continuous with `h(0) <= 0` (equality only in the unreachable
`κ -> ∞` limit) and `h(c_{fresh}) >= 0` (equality only at `κ = 0`, no
reaction — handled directly, without calling the root finder, since `h`
is then identically zero and has no isolated bracket), so `(0, c_{fresh})`
always brackets exactly one root, found by bisection.

Returns a named tuple `(c_in, c_out, p)`, with `p = 1 - c_out/c_fresh` the
overall extent of reaction against the fresh feed.
"""
function pfr_step_growth_self_catalyzed_recycle(τ, R, k, c_fresh)
    τ_hold = τ / (1 + R)
    κ = k * τ_hold
    if κ == 0
        return (c_in=c_fresh, c_out=c_fresh, p=0.0)
    end
    c_in_of(x) = (c_fresh + R * x) / (1 + R)
    h(x) = x - c_in_of(x) / sqrt(1 + 2 * κ * c_in_of(x)^2)
    c_out = find_zero(h, (0.0, c_fresh), Bisection())
    c_in = c_in_of(c_out)
    return (c_in=c_in, c_out=c_out, p=1 - c_out / c_fresh)
end
