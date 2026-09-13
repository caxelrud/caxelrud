"""
A small, general-purpose flowsheet framework: streams, mixers, splitters,
and the existing reactor unit operations wired together into arbitrary
topologies, including ones with recycle loops.

`reactor_recycle.jl` solved *one specific* topology (one PFR, one recycle
loop around it) by hand, deriving a closed form for each of the four
kinetic schemes. That doesn't generalize: two reactors with independent
recycle loops, a splitter feeding two parallel reactors that remix
downstream, or a reactor train with recycle around just one stage, all
need their own hand derivation under that approach. This module instead
provides the general *mechanism* standard process simulators use
(including commercial tools like Aspen): represent the flowsheet as
ordinary function composition over [`Stream`](@ref) values (mixers and
splitters are provided; reactors are the existing `cstr_*`/`pfr_*`
functions, extended here to accept and return `Stream`s via new methods),
and close any recycle loop by *tear-stream* convergence — guess the torn
stream, propagate it around the loop, and iterate the guess toward its own
image under that map ([`solve_tear`](@ref)). This is the standard
"sequential-modular with tear streams" method taught in any process design
course (e.g. Seader, Henley & Roper, *Separation Process Principles*) —
not something specific to this package's own derivations, unlike the
Sanchez-Lacombe mixture chemical-potential question this package
deliberately stays out of.

Because it's the same method under the hood, the general solver is
checked directly against `reactor_recycle.jl`'s hand-derived closed forms
in the test suite: assembling the *same* single-PFR-with-recycle topology
out of `mix`/`split_stream`/`pfr_free_radical` and solving it with
[`solve_tear`](@ref) must reproduce [`pfr_free_radical_recycle`](@ref) (and
the analogous coordination/step-growth functions) to numerical precision.
"""

"""
    Stream(Q, x)

A process stream: `Q` is the volumetric flow rate (any consistent unit —
the reactor unit operations only ever use flow-rate *ratios*), and `x` is
a `NamedTuple` of species concentrations, using the same field names the
corresponding scalar reactor functions use (`I`/`M` for free-radical,
`C_star`/`M` for coordination, `c` for step-growth). Constant density is
assumed throughout (as it already is in every `cstr_*`/`pfr_*` function,
whose `τ = V/Q` presumes a single flow rate in and out) — a reactor does
not change `Q`, only `x`; a splitter changes `Q` but not `x`; a mixer
combines both.
"""
struct Stream{T<:NamedTuple}
    Q::Float64
    x::T
end
Stream(Q, x::NamedTuple) = Stream(Float64(Q), x)

"""
    mix(streams::Stream...)

Combine two or more streams into one, by flow-weighted averaging of every
concentration field (all input streams must share the same field names):

``Q = \\sum_i Q_i, \\qquad x_j = \\dfrac{\\sum_i Q_i x_{i,j}}{\\sum_i Q_i}``

This is the physical statement that mixing conserves both total volumetric
flow (constant density) and species mass/moles.
"""
function mix(streams::Stream...)
    length(streams) >= 1 || throw(ArgumentError("mix needs at least one stream"))
    Qtot = sum(s.Q for s in streams)
    Qtot > 0 || throw(ArgumentError("total flow rate must be positive"))
    x = map(k -> sum(s.Q * getfield(s.x, k) for s in streams) / Qtot, fieldnames(typeof(streams[1].x)))
    return Stream(Qtot, NamedTuple{fieldnames(typeof(streams[1].x))}(x))
end

"""
    split_stream(stream::Stream, fraction)

Split `stream` into two streams of the *same* composition (a splitter
cannot change composition, only flow), with `fraction` of the flow going
to the first output and `1 - fraction` to the second:

``(Q_1, x) = (fraction \\cdot Q,\\ x), \\qquad (Q_2, x) = ((1-fraction)\\cdot Q,\\ x)``

`0 <= fraction <= 1` is required. Returns the pair `(stream1, stream2)`.

Named `split_stream` rather than `split` to avoid shadowing `Base.split`
(string splitting) — `using PolyRigorous` would otherwise make plain
`split` ambiguous in any script that also uses `Base.split`.
"""
function split_stream(stream::Stream, fraction)
    0 <= fraction <= 1 || throw(ArgumentError("fraction must be in [0, 1]"))
    return (Stream(fraction * stream.Q, stream.x), Stream((1 - fraction) * stream.Q, stream.x))
end

_field_close(a, b; atol, rtol) = abs(a - b) <= atol + rtol * max(abs(a), abs(b))

function _stream_close(s1::Stream, s2::Stream; atol, rtol)
    _field_close(s1.Q, s2.Q; atol=atol, rtol=rtol) || return false
    return all(_field_close(getfield(s1.x, k), getfield(s2.x, k); atol=atol, rtol=rtol) for k in fieldnames(typeof(s1.x)))
end

function _stream_step(s_old::Stream, s_new::Stream, damping)
    Q = s_old.Q + damping * (s_new.Q - s_old.Q)
    x = map((a, b) -> a + damping * (b - a), s_old.x, s_new.x)
    return Stream(Q, x)
end

"""
    solve_tear(f, x0::Stream; atol=1e-12, rtol=1e-9, maxiter=10_000, damping=1.0)

Converge a flowsheet's tear stream by successive substitution. `f` is the
"go around the loop once" map: given a guess for the torn stream, mix it
with whatever fresh feed(s) it recombines with, push the result through
the unit operation(s) that close the loop, use [`split_stream`](@ref) to
carve off the new guess for the same stream, and return it. The fixed point `x = f(x)` is exactly the
condition the real flowsheet satisfies at steady state.

Iterates `x_{n+1} = x_n + damping \\cdot (f(x_n) - x_n)` — plain
successive substitution at `damping = 1` (the default; this is exactly
what commercial process simulators try first), or damped/under-relaxed
for `damping < 1` if a particular loop doesn't contract on its own.
Convergence is judged on *every* field of the stream (flow rate and each
concentration) independently, combining an absolute and relative
tolerance (`abs(new - old) <= atol + rtol * max(abs(new), abs(old))`) so
it stays meaningful for a field that happens to sit near zero.

Returns a named tuple `(stream, converged, iterations)`. `converged` is
checked directly rather than assumed — a loop that doesn't contract under
plain successive substitution needs `damping < 1`, exactly as it would in
a real flowsheeting tool.
"""
function solve_tear(f, x0::Stream; atol=1e-12, rtol=1e-9, maxiter=10_000, damping=1.0)
    x = x0
    for i in 1:maxiter
        x_new = f(x)
        if _stream_close(x, x_new; atol=atol, rtol=rtol)
            return (stream=x_new, converged=true, iterations=i)
        end
        x = _stream_step(x, x_new, damping)
    end
    return (stream=x, converged=false, iterations=maxiter)
end

# ----------------------------------------------------------------------------
# Stream-aware methods for the existing reactor unit operations. Each wraps
# the already-tested scalar function(s) it's named after; none of these
# introduce new physics.
# ----------------------------------------------------------------------------

function cstr_free_radical(τ, k_p, k_d, k_t, f, s_in::Stream)
    res = cstr_free_radical(τ, k_p, k_d, k_t, f, s_in.x.I, s_in.x.M)
    return Stream(s_in.Q, (I=res.I, M=res.M))
end

"""
    pfr_free_radical(τ, k_p, k_d, k_t, f, s_in::Stream)

`Stream`-based PFR for free-radical polymerization: outlet monomer via
[`monomer_concentration`](@ref) as usual, and outlet initiator from its own
batch decay `I_out = I_in exp(-k_d τ)` (the same closed form used inside
[`pfr_free_radical_recycle`](@ref)) — needed here, unlike the scalar
`pfr_free_radical`, because a `Stream`'s outlet composition may need to
feed a *further* downstream unit operation.
"""
function pfr_free_radical(τ, k_p, k_d, k_t, f, s_in::Stream)
    I_out = s_in.x.I * exp(-k_d * τ)
    M_out = monomer_concentration(τ, k_p, k_d, k_t, f, s_in.x.I, s_in.x.M)
    return Stream(s_in.Q, (I=I_out, M=M_out))
end

function cstr_coordination(τ, k_p, k_t, s_in::Stream)
    res = cstr_coordination(τ, k_p, k_t, s_in.x.C_star, s_in.x.M)
    return Stream(s_in.Q, (C_star=res.C_star, M=res.M))
end

"""
    pfr_coordination(τ, k_p, k_t, s_in::Stream)

`Stream`-based PFR for coordination polymerization, analogous to
[`pfr_free_radical(τ, k_p, k_d, k_t, f, s_in::Stream)`](@ref): outlet
active-site concentration from its own batch decay
`C*_out = C*_in exp(-k_t τ)`, outlet monomer via
[`monomer_concentration_coordination`](@ref).
"""
function pfr_coordination(τ, k_p, k_t, s_in::Stream)
    C_star_out = s_in.x.C_star * exp(-k_t * τ)
    M_out = monomer_concentration_coordination(τ, k_p, s_in.x.C_star, k_t, s_in.x.M)
    return Stream(s_in.Q, (C_star=C_star_out, M=M_out))
end

function cstr_step_growth_external_catalyst(τ, k, s_in::Stream)
    res = cstr_step_growth_external_catalyst(τ, k, s_in.x.c)
    return Stream(s_in.Q, (c=res.c,))
end

function cstr_step_growth_self_catalyzed(τ, k, s_in::Stream)
    res = cstr_step_growth_self_catalyzed(τ, k, s_in.x.c)
    return Stream(s_in.Q, (c=res.c,))
end

function pfr_step_growth_external_catalyst(τ, k, s_in::Stream)
    p = extent_reaction_external_catalyst(k, s_in.x.c, τ)
    return Stream(s_in.Q, (c=s_in.x.c * (1 - p),))
end

function pfr_step_growth_self_catalyzed(τ, k, s_in::Stream)
    p = extent_reaction_self_catalyzed(k, s_in.x.c, τ)
    return Stream(s_in.Q, (c=s_in.x.c * (1 - p),))
end
