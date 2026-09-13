"""
Multi-zone reactors with staged initiator injection: the two classic
high-pressure (150-300 MPa) LDPE free-radical process configurations not
covered by `reactors.jl`'s plain single-pass or trains — a **high-pressure
autoclave** (one or more well-mixed zones, industrially with fresh
initiator injected at each zone to sustain the radical population as the
previous zone's charge decays) and a **tubular reactor** (a single long
PFR with multiple initiator injection points along its length, for the
same reason).

Convention: `τs` is a vector of per-stage/per-segment residence times
(zones for the autoclave, segments between injection points for the
tubular reactor), and `ΔI` a same-length vector of *additional* initiator
concentration injected at the start of each stage — on top of whatever
initiator concentration survived (decayed) from the previous stage, not
replacing it. Monomer is not re-dosed (the injection streams are treated
as small side streams that add initiator without materially diluting the
main monomer flow — the standard simplifying assumption for this kind of
process model). Setting `ΔI = [I0, 0, 0, ...]` (all the initiator charged
up front, nothing thereafter) recovers a plain (non-staged) CSTR train or
PFR exactly — checked in the test suite, along with the property that
splitting a plain PFR into extra zero-injection segments never changes
the result (confirming the staging machinery itself introduces no error,
independent of the injection feature).
"""

"""
    cstr_train_free_radical_staged(τs, k_p, k_d, k_t, f, ΔI, M_in)

Multi-zone high-pressure autoclave: `n = length(τs)` well-mixed zones in
series (each an ideal CSTR, [`cstr_free_radical`](@ref)), with fresh
initiator `ΔI[i]` injected into zone `i`'s own inlet on top of the
initiator concentration carried over from zone `i-1`'s outlet (zero for
the first zone).

Returns a named tuple `(I, M, conversion)` for the last zone's outlet,
`conversion` against the overall feed `M_in`.
"""
function cstr_train_free_radical_staged(τs, k_p, k_d, k_t, f, ΔI, M_in)
    length(ΔI) == length(τs) || throw(ArgumentError("ΔI and τs must have the same length"))
    I, M = 0.0, M_in
    for (τ, δI) in zip(τs, ΔI)
        I_in = I + δI
        res = cstr_free_radical(τ, k_p, k_d, k_t, f, I_in, M)
        I, M = res.I, res.M
    end
    return (I=I, M=M, conversion=1 - M / M_in)
end

"""
    pfr_train_free_radical_staged(τs, k_p, k_d, k_t, f, ΔI, M_in)

Tubular reactor with staged initiator injection: `n = length(τs)`
PFR segments in series (each kinetically a batch reactor run for time
`τs[i]`, as in [`pfr_free_radical`](@ref)), with fresh initiator `ΔI[i]`
injected at the *start* of segment `i` on top of the initiator
concentration surviving (via batch decay) from segment `i-1`'s outlet.

Returns a named tuple `(I, M, conversion)` for the last segment's outlet,
`conversion` against the overall feed `M_in`.
"""
function pfr_train_free_radical_staged(τs, k_p, k_d, k_t, f, ΔI, M_in)
    length(ΔI) == length(τs) || throw(ArgumentError("ΔI and τs must have the same length"))
    I, M = 0.0, M_in
    for (τ, δI) in zip(τs, ΔI)
        I_in = I + δI
        I = I_in * exp(-k_d * τ)
        M = monomer_concentration(τ, k_p, k_d, k_t, f, I_in, M)
    end
    return (I=I, M=M, conversion=1 - M / M_in)
end

"""
    ideal_gas_concentration(P, T)

Ideal-gas monomer concentration `[M] = P/(RT)` (mol/L) for a monomer
partial pressure `P` (MPa) at temperature `T` (K) — the standard way to
convert a gas-phase reactor's monomer partial pressure into the
concentration this package's kinetics functions expect (`k_p [M] [C*]`,
etc.), e.g. for modeling gas-phase fluidized-bed polyolefin reactors with
[`cstr_coordination`](@ref)/[`pfr_coordination`](@ref).
"""
function ideal_gas_concentration(P, T)
    R_GAS = 8.314462618  # J/(mol K)
    return (P * 1e6) / (R_GAS * T) / 1000  # Pa -> mol/m^3 -> mol/L
end
