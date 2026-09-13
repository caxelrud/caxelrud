"""
Coordination (Ziegler-Natta / metallocene-type) polymerization kinetics.

Convention: `C_star` is the concentration of active (chain-carrying)
catalyst sites, `M` the monomer concentration, `k_p` the propagation rate
constant.

Unlike free-radical termination — bimolecular, and it consumes two
growing chains at once — coordination polymerization is "living-like":
each chain-transfer event releases exactly one dead polymer molecule and
*regenerates* a new growing chain at the same site (so the active-site
count doesn't change), and each true-termination event releases one dead
polymer molecule and permanently kills the site. Either way, exactly one
dead chain is created per release event.

This lets `Xn` be *derived* from conservation of chain count rather than
needing an uncertain closed-form textbook formula: over any period,

``X_n = \\dfrac{\\text{monomer units consumed by propagation}}{\\text{dead polymer molecules released}}
       = \\dfrac{R_p}{\\text{rate of all chain-releasing events}}``

— see [`Xn_coordination`](@ref).
"""

"""
    coordination_propagation_rate(k_p, C_star, M)

Rate of monomer consumption by propagation at active coordination sites,
`Rp = k_p [C*] [M]` — the coordination-polymerization analogue of
[`propagation_rate`](@ref).
"""
coordination_propagation_rate(k_p, C_star, M) = k_p * C_star * M

"""
    transfer_rate_constant(; k_trM=0.0, M=0.0, k_trH=0.0, H2=0.0, k_tr0=0.0, k_t=0.0)

Convenience helper building the pseudo-first-order total chain-release
rate constant for [`Xn_coordination`](@ref) from the individual named
mechanisms commonly distinguished in coordination polymerization
(Ziegler-Natta, metallocene):

- `k_trM * M`: chain transfer to monomer.
- `k_trH * H2`: chain transfer to hydrogen (used industrially as a
  molecular-weight control agent).
- `k_tr0`: spontaneous transfer (e.g. β-hydride elimination), independent
  of any other species' concentration.
- `k_t`: true (irreversible) catalyst site termination.

Every term defaults to zero, so pass only the mechanisms relevant to your
system. Note that only `k_t` reduces the active-site count over time (see
[`monomer_concentration_coordination`](@ref)); the transfer mechanisms
release a dead chain but immediately start a new one at the same site.
"""
function transfer_rate_constant(; k_trM=0.0, M=0.0, k_trH=0.0, H2=0.0, k_tr0=0.0, k_t=0.0)
    return k_trM * M + k_trH * H2 + k_tr0 + k_t
end

"""
    Xn_coordination(k_p, M, k_release_total)

Number-average degree of polymerization for coordination (living-like)
chain-growth polymerization:

``X_n = \\dfrac{k_p [M]}{k_{release,total}}``

`k_release_total` is the pseudo-first-order total rate constant for all
chain-releasing pathways combined (see [`transfer_rate_constant`](@ref)
to build it up from individual mechanisms) — the active-site concentration
`[C*]` cancels out of both the propagation rate and the release rate, so
it doesn't appear here. Like [`kinetic_chain_length`](@ref) for
free-radical polymerization, this is an *instantaneous* quantity: the
average length of chains being formed right now, at these conditions —
not a cumulative average over a batch run in which conditions change.
"""
function Xn_coordination(k_p, M, k_release_total)
    k_release_total > 0 || throw(ArgumentError("k_release_total must be positive"))
    return k_p * M / k_release_total
end

"""
    monomer_concentration_coordination(t, k_p, C0_star, k_t, M0)

Monomer concentration at time `t` in an isothermal batch reactor running
coordination polymerization, where the active-site concentration decays
by true termination alone, `[C*](t) = C0_star exp(-k_t t)` (chain-transfer
events do *not* change `[C*]`, by construction — see the module docs), so

``\\frac{d[M]}{dt} = -k_p [C^*](t) [M] = -k_p C^*_0 e^{-k_t t} [M]``

integrates in closed form to

``\\ln\\frac{[M]_0}{[M](t)} = \\frac{k_p C^*_0}{k_t}\\bigl(1 - e^{-k_t t}\\bigr)``

with the `k_t -> 0` (non-deactivating, fully "living") limit handled
separately below as `k_p C^*_0 t`, since the formula above is `0/0` there
(this is the same kind of numerical pitfall as the free-radical batch
solution's `k_d -> 0` limit — caught by direct testing, not assumed away).
"""
function monomer_concentration_coordination(t, k_p, C0_star, k_t, M0)
    ln_ratio = if k_t == 0
        k_p * C0_star * t
    else
        (k_p * C0_star / k_t) * (1 - exp(-k_t * t))
    end
    return M0 * exp(-ln_ratio)
end

"""
    conversion_coordination(t, k_p, C0_star, k_t, M0)

Fractional monomer conversion `1 - [M](t)/[M]_0` at time `t`, from
[`monomer_concentration_coordination`](@ref).
"""
conversion_coordination(t, k_p, C0_star, k_t, M0) = 1 - monomer_concentration_coordination(t, k_p, C0_star, k_t, M0) / M0
