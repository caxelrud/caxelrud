using Test
using PolyRigorous

@testset "PolyRigorous" begin

@testset "components" begin
    ps = species("polystyrene")
    tol = species("toluene")
    @test ps.kind == :polymer
    @test tol.kind == :solvent

    N = degree_of_polymerization(ps, tol)
    @test N > 1  # a 100 kg/mol chain is much bigger than a solvent molecule

    r_tol = segment_number(tol)
    @test r_tol > 0

    @test_throws KeyError species("unobtainium")

    # Polyethylene (both grades) and polypropylene must be present and
    # give sensible, solvable Sanchez-Lacombe PVT behavior.
    ldpe = species("polyethylene_ldpe")
    hdpe = species("polyethylene")
    pp = species("polypropylene")
    for sp in (ldpe, hdpe, pp)
        @test sp.kind == :polymer
        @test segment_number(sp) > 0
        ρ = density(298.15, 0.1, sp)
        @test 0 < ρ < sp.rhostar
    end
end

@testset "flory-huggins: limits and known values" begin
    # Athermal, equal-size (N = 1) binary: ln(a1) reduces to ideal-solution
    # (Raoult's-law-like) behavior ln(a1) = ln(φ1).
    φ₂ = 0.3
    @test isapprox(PolyRigorous.ln_activity_solvent(φ₂, 1.0, 0.0), log(1 - φ₂); atol=1e-12)

    # As N -> infinity at fixed φ₂, χ (the (1-1/N) term saturates)
    ln_a1_bigN = PolyRigorous.ln_activity_solvent(φ₂, 1e8, 0.5)
    ln_a1_inftyN = log(1 - φ₂) + φ₂ + 0.5 * φ₂^2
    @test isapprox(ln_a1_bigN, ln_a1_inftyN; atol=1e-6)

    # Symmetry check for N=1: solvent and polymer are the same size, so
    # relabeling which species is "solvent" (φ₂ <-> 1-φ₂) should map
    # ln_activity_solvent onto ln_activity_polymer.
    @test isapprox(
        PolyRigorous.ln_activity_solvent(φ₂, 1.0, 0.7),
        PolyRigorous.ln_activity_polymer(1 - φ₂, 1.0, 0.7);
        atol=1e-12,
    )
end

@testset "flory-huggins: critical point and spinodal" begin
    N = 1.0
    crit = critical_point(N)
    @test isapprox(crit.φ₂_c, 0.5; atol=1e-12)
    @test isapprox(crit.χ_c, 2.0; atol=1e-12)

    N2 = 1000.0
    crit2 = critical_point(N2)
    @test isapprox(crit2.φ₂_c, 1 / (1 + sqrt(N2)); atol=1e-12)

    # The spinodal touches the binodal exactly at the critical point:
    # chi_spinodal(φ₂_c) == χ_c.
    @test isapprox(chi_spinodal(crit2.φ₂_c, N2), crit2.χ_c; atol=1e-9)
end

@testset "phase equilibrium: binodal brackets the critical point" begin
    N = 50.0
    crit = critical_point(N)
    χ = 1.05 * crit.χ_c
    res = binodal_pair(N, χ)
    @test res.converged
    @test res.φ₂a < crit.φ₂_c < res.φ₂b

    # Both branches should satisfy the equal-activity conditions to high
    # precision.
    a1a = PolyRigorous.ln_activity_solvent(res.φ₂a, N, χ)
    a1b = PolyRigorous.ln_activity_solvent(res.φ₂b, N, χ)
    @test isapprox(a1a, a1b; atol=1e-6)

    a2a = PolyRigorous.ln_activity_polymer(res.φ₂a, N, χ)
    a2b = PolyRigorous.ln_activity_polymer(res.φ₂b, N, χ)
    @test isapprox(a2a, a2b; atol=1e-6)
end

@testset "phase equilibrium: binodal_curve and spinodal_curve run and stay ordered" begin
    N = 80.0
    curve = binodal_curve(N; npoints=15)
    @test length(curve) > 5
    for pt in curve
        @test pt.φ₂a < pt.φ₂b
    end

    spin = spinodal_curve(N; npoints=50)
    @test length(spin) == 50
    @test all(pt.χ > 0 for pt in spin)
end

@testset "sanchez-lacombe: EOS residual and PVT sanity" begin
    tol = species("toluene")
    T, P = 298.15, 0.1  # K, MPa (~1 atm)

    ρ̃ = PolyRigorous.reduced_density(T, P, tol)
    @test 0 < ρ̃ < 1

    r = segment_number(tol)
    resid = sl_eos_residual(ρ̃, T / tol.Tstar, P / tol.Pstar, r)
    @test isapprox(resid, 0.0; atol=1e-8)

    ρ = density(T, P, tol)
    @test 0 < ρ < tol.rhostar

    # Increasing pressure at fixed T should compress the fluid (density
    # goes up).
    ρ_hi = density(T, 50.0, tol)
    @test ρ_hi > ρ

    # Increasing temperature at fixed P should expand the fluid (density
    # goes down) below the polymer's/solvent's spinodal-like limits.
    ρ_hot = density(T + 20, P, tol)
    @test ρ_hot < ρ

    α = thermal_expansion_coefficient(T, P, tol)
    @test α > 0

    κ = isothermal_compressibility(T, P, tol)
    @test κ > 0
end

@testset "sanchez-lacombe: binary mixture PVT" begin
    tol = species("toluene")
    ps = species("polystyrene")

    @test_throws ArgumentError sl_mixing_rules(tol, 0.6, ps, 0.6)

    # Pure-component limits: at w1=1 (resp. w1=0) the mixing rules must
    # reduce exactly to component 1's (resp. component 2's) own
    # parameters, since phi1=1, phi2=0 (resp. the reverse) makes every
    # cross term drop out.
    mix1 = sl_mixing_rules(tol, 1.0, ps, 0.0)
    @test isapprox(mix1.Tstar, tol.Tstar; rtol=1e-10)
    @test isapprox(mix1.Pstar, tol.Pstar; rtol=1e-10)
    @test isapprox(mix1.rhostar, tol.rhostar; rtol=1e-10)
    @test isapprox(mix1.M, tol.M; rtol=1e-10)
    @test isapprox(mix1.r, segment_number(tol); rtol=1e-10)

    mix2 = sl_mixing_rules(tol, 0.0, ps, 1.0)
    @test isapprox(mix2.Tstar, ps.Tstar; rtol=1e-10)
    @test isapprox(mix2.rhostar, ps.rhostar; rtol=1e-10)
    @test isapprox(mix2.r, segment_number(ps); rtol=1e-10)

    # M_mix is constructed specifically so that segment_number, applied to
    # the resulting virtual Species, reproduces the mixing rule's own r --
    # this is the algebraic identity sl_mixture_species relies on to reuse
    # the pure-component EOS machinery unmodified.
    w1 = 0.3
    mix = sl_mixing_rules(tol, w1, ps, 1 - w1)
    sp_mix = sl_mixture_species(tol, w1, ps, 1 - w1)
    @test isapprox(segment_number(sp_mix), mix.r; rtol=1e-10)
    @test isapprox(sp_mix.Tstar, mix.Tstar; rtol=1e-12)
    @test isapprox(sp_mix.Pstar, mix.Pstar; rtol=1e-12)
    @test isapprox(sp_mix.rhostar, mix.rhostar; rtol=1e-12)

    # Consequently, PVT functions applied to the virtual mixture Species
    # must reduce to the pure-component result at the composition limits.
    T, P = 298.15, 0.1
    sp_pure1 = sl_mixture_species(tol, 1.0, ps, 0.0)
    @test isapprox(density(T, P, sp_pure1), density(T, P, tol); rtol=1e-8)
    sp_pure2 = sl_mixture_species(tol, 0.0, ps, 1.0)
    @test isapprox(density(T, P, sp_pure2), density(T, P, ps); rtol=1e-8)

    # At an intermediate composition, PVT functions should run cleanly and
    # give a sensible (finite, positive, sub-close-packed) result.
    ρ_mix = density(T, P, sp_mix)
    @test 0 < ρ_mix < sp_mix.rhostar

    # The binary interaction parameter k12 should have a real (nonzero)
    # effect on the mixed characteristic pressure for a pair with
    # different characteristic energies.
    mix_k0 = sl_mixing_rules(tol, w1, ps, 1 - w1; k12=0.0)
    mix_k = sl_mixing_rules(tol, w1, ps, 1 - w1; k12=0.05)
    @test !isapprox(mix_k0.Pstar, mix_k.Pstar; rtol=1e-6)
end

@testset "free-radical kinetics" begin
    k_p, k_d, k_t, f = 1e3, 1e-5, 1e7, 0.5
    I0, M0 = 0.01, 5.0

    # No reaction has happened yet at t=0.
    @test monomer_concentration(0.0, k_p, k_d, k_t, f, I0, M0) == M0
    @test conversion(0.0, k_p, k_d, k_t, f, I0, M0) == 0.0

    # Conversion increases monotonically in t and never exceeds 1 -- and, in
    # this dead-end (finite initiator charge) batch model, never reaches it
    # either: once the initiator is exhausted, propagation stops with some
    # monomer left unreacted. The limiting conversion as t -> infinity has a
    # known closed form (exp(-k_d t/2) -> 0 in the conversion formula).
    ts = [0.0, 100.0, 1e4, 1e6, 1e9]
    convs = [conversion(t, k_p, k_d, k_t, f, I0, M0) for t in ts]
    @test issorted(convs)
    @test all(0 .<= convs .< 1)
    conv_limit = 1 - exp(-(2 * k_p / k_d) * sqrt(f * k_d * I0 / k_t))
    @test isapprox(convs[end], conv_limit; atol=1e-9)

    # Two independent routes to the kinetic chain length should agree:
    # ν = Rp/Ri directly, vs. the closed-form kinetic_chain_length formula.
    Mrad = radical_concentration(f, k_d, I0, k_t)
    Ri = initiation_rate(f, k_d, I0)
    Rp = propagation_rate(k_p, M0, Mrad)
    ν_direct = Rp / Ri
    ν_formula = kinetic_chain_length(k_p, M0, f, k_d, I0, k_t)
    @test isapprox(ν_direct, ν_formula; rtol=1e-10)

    # Xn_mixed reduces to the pure combination/disproportionation cases at
    # its endpoints.
    ν = ν_formula
    @test isapprox(Xn_mixed(ν, 0.0), Xn_combination(ν); atol=1e-10)
    @test isapprox(Xn_mixed(ν, 1.0), Xn_disproportionation(ν); atol=1e-10)
    @test Xn_combination(ν) > Xn_disproportionation(ν)  # 2ν > ν
end

@testset "step-growth kinetics and the Flory distribution" begin
    k, c0 = 0.5, 1.0

    @test extent_reaction_external_catalyst(k, c0, 0.0) == 0.0
    @test extent_reaction_self_catalyzed(k, c0, 0.0) == 0.0

    # Both extents of reaction increase monotonically toward 1 (complete
    # reaction) as t -> infinity, but never reach or exceed it.
    ts = [0.0, 1.0, 10.0, 1e4, 1e8]
    p_ext = [extent_reaction_external_catalyst(k, c0, t) for t in ts]
    p_self = [extent_reaction_self_catalyzed(k, c0, t) for t in ts]
    @test issorted(p_ext) && all(0 .<= p_ext .< 1)
    @test issorted(p_self) && all(0 .<= p_self .< 1)
    @test p_ext[end] > 0.999
    @test p_self[end] > 0.999

    # Carothers equation: no reaction means no polymerization (Xn = 1);
    # the stoichiometrically-balanced case (r=1) is just 1/(1-p).
    @test carothers_Xn(0.0) == 1.0
    p = 0.95
    @test isapprox(carothers_Xn(p), 1 / (1 - p); atol=1e-12)
    @test isapprox(carothers_Xn(p; r=1.0), 1 / (1 - p); atol=1e-12)
    # A stoichiometric imbalance (r < 1) caps Xn below the balanced case.
    @test carothers_Xn(p; r=0.98) < carothers_Xn(p; r=1.0)

    # Flory "most probable" distribution: mole fractions and weight
    # fractions are each known analytically to sum to 1 over all chain
    # lengths (geometric series identities); check this numerically with a
    # truncated but very long sum.
    xs = 1:200_000
    mole_sum = sum(flory_mole_fraction(x, p) for x in xs)
    weight_sum = sum(flory_weight_fraction(x, p) for x in xs)
    @test isapprox(mole_sum, 1.0; atol=1e-6)
    @test isapprox(weight_sum, 1.0; atol=1e-3)

    # PDI = Xw/Xn = 1+p, and approaches (but never reaches) 2 as p -> 1.
    @test isapprox(flory_Xw(p) / carothers_Xn(p), flory_PDI(p); atol=1e-10)
    @test isapprox(flory_PDI(0.0), 1.0; atol=1e-12)
    @test flory_PDI(0.9999) < 2.0
    @test flory_PDI(0.9999) > 1.999
end

@testset "reactors: free-radical CSTR and PFR" begin
    k_p, k_d, k_t, f = 1e3, 1e-5, 1e7, 0.5
    I_in, M_in = 0.01, 5.0
    τ = 3600.0  # s

    # As τ -> 0, nothing has reacted yet.
    res0 = cstr_free_radical(1e-8, k_p, k_d, k_t, f, I_in, M_in)
    @test isapprox(res0.I, I_in; rtol=1e-6)
    @test isapprox(res0.M, M_in; rtol=1e-6)
    @test isapprox(res0.conversion, 0.0; atol=1e-6)

    res = cstr_free_radical(τ, k_p, k_d, k_t, f, I_in, M_in)
    @test 0 < res.conversion < 1

    # The closed-form CSTR solution must satisfy the steady-state species
    # balances it was derived from: (Cin - C)/τ = consumption rate of C.
    @test isapprox((I_in - res.I) / τ, k_d * res.I; rtol=1e-10)
    Mrad = radical_concentration(f, k_d, res.I, k_t)
    @test isapprox((M_in - res.M) / τ, propagation_rate(k_p, res.M, Mrad); rtol=1e-10)

    # PFR is just the batch solution at t=τ -- check the wrapper is wired
    # correctly, and that a single CSTR is less efficient than a PFR at
    # the same residence time (true for these positive-order, concentration-
    # decreasing-rate kinetics: a CSTR runs at its low outlet concentration
    # the whole time, whereas a PFR/batch starts at the high feed
    # concentration where the rate is fastest).
    pfr_res = pfr_free_radical(τ, k_p, k_d, k_t, f, I_in, M_in)
    @test isapprox(pfr_res.M, monomer_concentration(τ, k_p, k_d, k_t, f, I_in, M_in); rtol=1e-12)
    @test res.conversion < pfr_res.conversion
end

@testset "reactors: step-growth CSTR and PFR" begin
    k, c_in = 0.5, 1.0
    τ = 20.0

    # As τ -> 0, nothing has reacted yet.
    res0_ext = cstr_step_growth_external_catalyst(1e-8, k, c_in)
    @test isapprox(res0_ext.c, c_in; rtol=1e-6)
    @test isapprox(res0_ext.p, 0.0; atol=1e-6)
    res0_self = cstr_step_growth_self_catalyzed(1e-8, k, c_in)
    @test isapprox(res0_self.c, c_in; rtol=1e-6)
    @test isapprox(res0_self.p, 0.0; atol=1e-6)

    res_ext = cstr_step_growth_external_catalyst(τ, k, c_in)
    @test 0 < res_ext.p < 1
    # Steady-state balance: (c_in - c)/τ = k c^2.
    @test isapprox((c_in - res_ext.c) / τ, k * res_ext.c^2; rtol=1e-10)

    res_self = cstr_step_growth_self_catalyzed(τ, k, c_in)
    @test 0 < res_self.p < 1
    # Steady-state balance: (c_in - c)/τ = k c^3.
    @test isapprox((c_in - res_self.c) / τ, k * res_self.c^3; rtol=1e-10)

    # PFR wrapper wiring, and the same CSTR-vs-PFR efficiency ordering as
    # the free-radical case.
    pfr_ext = pfr_step_growth_external_catalyst(τ, k, c_in)
    @test isapprox(pfr_ext.p, extent_reaction_external_catalyst(k, c_in, τ); rtol=1e-12)
    @test res_ext.p < pfr_ext.p

    pfr_self = pfr_step_growth_self_catalyzed(τ, k, c_in)
    @test isapprox(pfr_self.p, extent_reaction_self_catalyzed(k, c_in, τ); rtol=1e-12)
    @test res_self.p < pfr_self.p
end

@testset "reactors: CSTR trains converge to PFR as stage count grows" begin
    # A single-stage "train" must reduce exactly to a plain CSTR.
    k_p, k_d, k_t, f = 1e3, 1e-5, 1e7, 0.5
    I_in, M_in = 0.01, 5.0
    τ_total = 3600.0

    single = cstr_free_radical(τ_total, k_p, k_d, k_t, f, I_in, M_in)
    train1 = cstr_train_free_radical(1, τ_total, k_p, k_d, k_t, f, I_in, M_in)
    @test isapprox(train1.M, single.M; rtol=1e-12)

    # The vector-of-τs form must agree with the (n, τ) convenience form
    # when all stages are equal.
    n = 5
    train_vec = cstr_train_free_radical(fill(τ_total / n, n), k_p, k_d, k_t, f, I_in, M_in)
    train_n = cstr_train_free_radical(n, τ_total / n, k_p, k_d, k_t, f, I_in, M_in)
    @test isapprox(train_vec.M, train_n.M; rtol=1e-12)

    # Classic reactor-engineering result: at fixed *total* residence time,
    # a CSTR train's conversion increases monotonically with stage count
    # and converges to the PFR value as n -> infinity.
    pfr_res = pfr_free_radical(τ_total, k_p, k_d, k_t, f, I_in, M_in)
    ns = [1, 2, 5, 20, 100, 2000]
    convs = [cstr_train_free_radical(n, τ_total / n, k_p, k_d, k_t, f, I_in, M_in).conversion for n in ns]
    @test issorted(convs)
    @test all(convs .< pfr_res.conversion)
    @test isapprox(convs[end], pfr_res.conversion; atol=1e-3)

    # Same convergence property for step-growth, both kinetic orders.
    k, c_in = 0.5, 1.0
    pfr_ext = pfr_step_growth_external_catalyst(τ_total, k, c_in)
    convs_ext = [cstr_train_step_growth_external_catalyst(n, τ_total / n, k, c_in).p for n in ns]
    @test issorted(convs_ext)
    @test all(convs_ext .< pfr_ext.p)
    @test isapprox(convs_ext[end], pfr_ext.p; atol=1e-3)

    pfr_self = pfr_step_growth_self_catalyzed(τ_total, k, c_in)
    convs_self = [cstr_train_step_growth_self_catalyzed(n, τ_total / n, k, c_in).p for n in ns]
    @test issorted(convs_self)
    @test all(convs_self .< pfr_self.p)
    @test isapprox(convs_self[end], pfr_self.p; atol=1e-3)
end

@testset "coordination polymerization kinetics" begin
    k_p, M, C_star = 50.0, 2.0, 1e-4

    # transfer_rate_constant: defaults to zero, and each named mechanism
    # contributes additively.
    @test transfer_rate_constant() == 0.0
    k_trM, k_trH, k_tr0, k_t = 0.1, 0.05, 1e-3, 1e-4
    H2 = 0.01
    k_release = transfer_rate_constant(; k_trM=k_trM, M=M, k_trH=k_trH, H2=H2, k_tr0=k_tr0, k_t=k_t)
    @test isapprox(k_release, k_trM * M + k_trH * H2 + k_tr0 + k_t; atol=1e-14)

    @test_throws ArgumentError Xn_coordination(k_p, M, 0.0)

    Xn = Xn_coordination(k_p, M, k_release)
    @test Xn > 0

    # Self-consistency with the mass-balance derivation this formula comes
    # from: Xn * (rate of chain-releasing events) == propagation rate.
    Rp = coordination_propagation_rate(k_p, C_star, M)
    release_rate = k_release * C_star
    @test isapprox(Xn * release_rate, Rp; rtol=1e-12)

    # No termination at all (k_release_total -> 0) means infinitely long
    # chains -- consistent with the "living" limit.
    @test Xn_coordination(k_p, M, 1e-10) > Xn_coordination(k_p, M, 1.0)

    # Batch conversion: the non-deactivating (k_t=0, "living") and
    # deactivating (k_t>0) branches of monomer_concentration_coordination
    # must agree as k_t -> 0 (this is exactly the kind of 0/0 numerical
    # trap that broke the step-growth CSTR earlier -- checked directly,
    # not just assumed continuous).
    C0_star, M0, t = 1e-4, 5.0, 100.0
    M_living = monomer_concentration_coordination(t, k_p, C0_star, 0.0, M0)
    M_near_living = monomer_concentration_coordination(t, k_p, C0_star, 1e-10, M0)
    @test isapprox(M_living, M_near_living; rtol=1e-6)
    @test isapprox(M_living, M0 * exp(-k_p * C0_star * t); rtol=1e-12)

    @test conversion_coordination(0.0, k_p, C0_star, k_t, M0) == 0.0

    # Living (k_t=0): the active-site pool never decays, so given enough
    # time *all* the monomer is eventually consumed -- unlike free-radical
    # or a deactivating coordination catalyst, there's no "dead-end" limit.
    ts = [0.0, 1.0, 10.0, 100.0, 1e4]
    convs_living = [conversion_coordination(t, k_p, C0_star, 0.0, M0) for t in ts]
    @test issorted(convs_living)
    @test convs_living[end] > 0.999

    # Deactivating (k_t>0): conversion approaches a finite limit < 1, from
    # the same closed form used for the free-radical batch solution.
    convs_deactivating = [conversion_coordination(t, k_p, C0_star, k_t, M0) for t in ts]
    @test issorted(convs_deactivating)
    conv_limit = 1 - exp(-(k_p * C0_star / k_t))
    @test isapprox(convs_deactivating[end], conv_limit; atol=1e-6)
    @test convs_deactivating[end] < convs_living[end]
end

@testset "reactors: coordination CSTR, PFR, and train" begin
    k_p, k_t = 50.0, 1e-4
    C_star_in, M_in = 1e-4, 5.0
    τ = 100.0

    # As τ -> 0, nothing has reacted yet.
    res0 = cstr_coordination(1e-8, k_p, k_t, C_star_in, M_in)
    @test isapprox(res0.C_star, C_star_in; rtol=1e-6)
    @test isapprox(res0.M, M_in; rtol=1e-6)
    @test isapprox(res0.conversion, 0.0; atol=1e-6)

    res = cstr_coordination(τ, k_p, k_t, C_star_in, M_in)
    @test 0 < res.conversion < 1

    # Steady-state species balances: (Cin - C)/τ = consumption rate of C.
    @test isapprox((C_star_in - res.C_star) / τ, k_t * res.C_star; rtol=1e-10)
    @test isapprox((M_in - res.M) / τ, k_p * res.M * res.C_star; rtol=1e-10)

    # PFR is the batch solution at t=τ; CSTR is less efficient than PFR at
    # the same residence time (same ordering argument as free-radical).
    pfr_res = pfr_coordination(τ, k_p, C_star_in, k_t, M_in)
    @test isapprox(pfr_res.M, monomer_concentration_coordination(τ, k_p, C_star_in, k_t, M_in); rtol=1e-12)
    @test res.conversion < pfr_res.conversion

    # CSTR train: single stage reduces to a plain CSTR, and conversion
    # converges to the PFR value as stage count grows at fixed total τ.
    train1 = cstr_train_coordination(1, τ, k_p, k_t, C_star_in, M_in)
    @test isapprox(train1.M, res.M; rtol=1e-12)

    ns = [1, 2, 5, 20, 100, 2000]
    convs = [cstr_train_coordination(n, τ / n, k_p, k_t, C_star_in, M_in).conversion for n in ns]
    @test issorted(convs)
    @test all(convs .< pfr_res.conversion)
    @test isapprox(convs[end], pfr_res.conversion; atol=1e-3)
end

@testset "reactors: PFR with recycle" begin
    # --- Free-radical ---
    k_p, k_d, k_t, f = 1e3, 1e-5, 1e7, 0.5
    I_fresh, M_fresh = 0.01, 5.0
    τ = 3600.0

    # R=0 must reduce exactly to the plain PFR.
    rec0 = pfr_free_radical_recycle(τ, 0.0, k_p, k_d, k_t, f, I_fresh, M_fresh)
    plain = pfr_free_radical(τ, k_p, k_d, k_t, f, I_fresh, M_fresh)
    @test isapprox(rec0.M_out, plain.M, rtol=1e-10)
    @test isapprox(rec0.I_in, I_fresh; rtol=1e-12)
    @test isapprox(rec0.conversion, plain.conversion; rtol=1e-10)

    # Mass-balance self-consistency: the reported inlet concentrations must
    # satisfy the mixing equation, and the outlet must be the batch/PFR
    # evolution of the inlet over τ_hold = τ/(1+R).
    R = 2.5
    rec = pfr_free_radical_recycle(τ, R, k_p, k_d, k_t, f, I_fresh, M_fresh)
    τ_hold = τ / (1 + R)
    @test isapprox(rec.I_in, (I_fresh + R * rec.I_out) / (1 + R); rtol=1e-10)
    @test isapprox(rec.M_in, (M_fresh + R * rec.M_out) / (1 + R); rtol=1e-10)
    @test isapprox(rec.I_out, rec.I_in * exp(-k_d * τ_hold); rtol=1e-10)
    @test isapprox(rec.M_out, monomer_concentration(τ_hold, k_p, k_d, k_t, f, rec.I_in, rec.M_in); rtol=1e-8)

    # Large-R limit: a PFR with (near-)infinite recycle behaves like a CSTR
    # of the same nominal τ (classic reactor-engineering result).
    rec_bigR = pfr_free_radical_recycle(τ, 1e5, k_p, k_d, k_t, f, I_fresh, M_fresh)
    cstr_res = cstr_free_radical(τ, k_p, k_d, k_t, f, I_fresh, M_fresh)
    @test isapprox(rec_bigR.conversion, cstr_res.conversion; rtol=1e-3)

    # --- Coordination ---
    k_p2, k_t2 = 50.0, 1e-4
    C0_star_fresh, M_fresh2 = 1e-4, 5.0
    τ2 = 100.0

    rec0_coord = pfr_coordination_recycle(τ2, 0.0, k_p2, k_t2, C0_star_fresh, M_fresh2)
    plain_coord = pfr_coordination(τ2, k_p2, C0_star_fresh, k_t2, M_fresh2)
    @test isapprox(rec0_coord.M_out, plain_coord.M; rtol=1e-10)

    rec_bigR_coord = pfr_coordination_recycle(τ2, 1e5, k_p2, k_t2, C0_star_fresh, M_fresh2)
    cstr_res_coord = cstr_coordination(τ2, k_p2, k_t2, C0_star_fresh, M_fresh2)
    @test isapprox(rec_bigR_coord.conversion, cstr_res_coord.conversion; rtol=1e-3)

    # --- Step-growth, external catalyst ---
    k, c_fresh = 0.5, 1.0
    τ3 = 20.0

    rec0_ext = pfr_step_growth_external_catalyst_recycle(τ3, 0.0, k, c_fresh)
    plain_ext = pfr_step_growth_external_catalyst(τ3, k, c_fresh)
    @test isapprox(rec0_ext.p, plain_ext.p; rtol=1e-10)

    R2 = 3.0
    rec_ext = pfr_step_growth_external_catalyst_recycle(τ3, R2, k, c_fresh)
    τ3_hold = τ3 / (1 + R2)
    @test isapprox(rec_ext.c_in, (c_fresh + R2 * rec_ext.c_out) / (1 + R2); rtol=1e-10)
    p_ext_hold = extent_reaction_external_catalyst(k, rec_ext.c_in, τ3_hold)
    @test isapprox(rec_ext.c_out, rec_ext.c_in * (1 - p_ext_hold); rtol=1e-8)

    rec_bigR_ext = pfr_step_growth_external_catalyst_recycle(τ3, 1e5, k, c_fresh)
    cstr_res_ext = cstr_step_growth_external_catalyst(τ3, k, c_fresh)
    @test isapprox(rec_bigR_ext.p, cstr_res_ext.p; rtol=1e-3)

    # --- Step-growth, self-catalyzed ---
    rec0_self = pfr_step_growth_self_catalyzed_recycle(τ3, 0.0, k, c_fresh)
    plain_self = pfr_step_growth_self_catalyzed(τ3, k, c_fresh)
    @test isapprox(rec0_self.p, plain_self.p; rtol=1e-8)

    rec_self = pfr_step_growth_self_catalyzed_recycle(τ3, R2, k, c_fresh)
    @test isapprox(rec_self.c_in, (c_fresh + R2 * rec_self.c_out) / (1 + R2); rtol=1e-8)
    p_self_hold = extent_reaction_self_catalyzed(k, rec_self.c_in, τ3_hold)
    @test isapprox(rec_self.c_out, rec_self.c_in * (1 - p_self_hold); rtol=1e-8)

    rec_bigR_self = pfr_step_growth_self_catalyzed_recycle(τ3, 1e5, k, c_fresh)
    cstr_res_self = cstr_step_growth_self_catalyzed(τ3, k, c_fresh)
    @test isapprox(rec_bigR_self.p, cstr_res_self.p; rtol=1e-3)

    # No reaction (k=0): every recycle case must report zero conversion.
    @test pfr_free_radical_recycle(τ, R, 0.0, k_d, k_t, f, I_fresh, M_fresh).conversion == 0.0
    @test isapprox(pfr_coordination_recycle(τ2, R, 0.0, k_t2, C0_star_fresh, M_fresh2).conversion, 0.0; atol=1e-12)
    @test pfr_step_growth_external_catalyst_recycle(τ3, R2, 0.0, c_fresh).p == 0.0
    @test pfr_step_growth_self_catalyzed_recycle(τ3, R2, 0.0, c_fresh).p == 0.0
end

@testset "flowsheet: Stream, mix, split_stream" begin
    s1 = Stream(2.0, (M=1.0,))
    s2 = Stream(3.0, (M=4.0,))
    m = mix(s1, s2)
    @test m.Q == 5.0
    @test isapprox(m.x.M, (2 * 1.0 + 3 * 4.0) / 5.0; atol=1e-12)

    # mix with more than two streams: weights must still sum correctly.
    s3 = Stream(5.0, (M=10.0,))
    m3 = mix(s1, s2, s3)
    @test m3.Q == 10.0
    @test isapprox(m3.x.M, (2 * 1.0 + 3 * 4.0 + 5 * 10.0) / 10.0; atol=1e-12)

    a, b = split_stream(Stream(10.0, (M=7.0,)), 0.3)
    @test isapprox(a.Q, 3.0; atol=1e-12) && isapprox(b.Q, 7.0; atol=1e-12)
    @test a.x.M == 7.0 && b.x.M == 7.0  # a splitter never changes composition
    @test_throws ArgumentError split_stream(Stream(10.0, (M=7.0,)), 1.5)
end

@testset "flowsheet: general solve_tear reproduces reactor_recycle closed forms" begin
    # For each of the four kinetic schemes, assemble the *same* single-PFR-
    # with-recycle topology out of mix/split_stream/pfr_* and solve it with
    # the general tear-stream solver -- it must reproduce the hand-derived
    # closed form in reactor_recycle.jl to numerical precision. This ties
    # the general mechanism directly to results already verified above.
    Q_fresh = 1.0
    R = 2.5

    # --- Free-radical ---
    k_p, k_d, k_t, f = 1e3, 1e-5, 1e7, 0.5
    I_fresh, M_fresh = 0.01, 5.0
    τ = 3600.0
    τ_hold = τ / (1 + R)
    s_fresh = Stream(Q_fresh, (I=I_fresh, M=M_fresh))

    function loop_frk(recycle)
        mixed = mix(s_fresh, recycle)
        out = pfr_free_radical(τ_hold, k_p, k_d, k_t, f, mixed)
        new_recycle, _ = split_stream(out, R / (1 + R))
        return new_recycle
    end
    res_frk = solve_tear(loop_frk, Stream(R * Q_fresh, s_fresh.x))
    @test res_frk.converged

    mixed_final = mix(s_fresh, res_frk.stream)
    out_final = pfr_free_radical(τ_hold, k_p, k_d, k_t, f, mixed_final)
    _, product_final = split_stream(out_final, R / (1 + R))

    closed_frk = pfr_free_radical_recycle(τ, R, k_p, k_d, k_t, f, I_fresh, M_fresh)
    @test isapprox(product_final.x.M, closed_frk.M_out; rtol=1e-6)
    @test isapprox(mixed_final.x.I, closed_frk.I_in; rtol=1e-6)
    @test isapprox(out_final.x.I, closed_frk.I_out; rtol=1e-6)

    # --- Coordination ---
    k_p2, k_t2 = 50.0, 1e-4
    C0_star_fresh, M_fresh2 = 1e-4, 5.0
    τ2 = 100.0
    τ2_hold = τ2 / (1 + R)
    s_fresh2 = Stream(Q_fresh, (C_star=C0_star_fresh, M=M_fresh2))

    function loop_coord(recycle)
        mixed = mix(s_fresh2, recycle)
        out = pfr_coordination(τ2_hold, k_p2, k_t2, mixed)
        new_recycle, _ = split_stream(out, R / (1 + R))
        return new_recycle
    end
    res_coord = solve_tear(loop_coord, Stream(R * Q_fresh, s_fresh2.x))
    @test res_coord.converged

    mixed2_final = mix(s_fresh2, res_coord.stream)
    out2_final = pfr_coordination(τ2_hold, k_p2, k_t2, mixed2_final)
    _, product2_final = split_stream(out2_final, R / (1 + R))

    closed_coord = pfr_coordination_recycle(τ2, R, k_p2, k_t2, C0_star_fresh, M_fresh2)
    @test isapprox(product2_final.x.M, closed_coord.M_out; rtol=1e-6)

    # --- Step-growth, external catalyst ---
    k, c_fresh = 0.5, 1.0
    τ3 = 20.0
    τ3_hold = τ3 / (1 + R)
    s_fresh3 = Stream(Q_fresh, (c=c_fresh,))

    function loop_ext(recycle)
        mixed = mix(s_fresh3, recycle)
        out = pfr_step_growth_external_catalyst(τ3_hold, k, mixed)
        new_recycle, _ = split_stream(out, R / (1 + R))
        return new_recycle
    end
    res_ext = solve_tear(loop_ext, Stream(R * Q_fresh, s_fresh3.x))
    @test res_ext.converged

    mixed3_final = mix(s_fresh3, res_ext.stream)
    out3_final = pfr_step_growth_external_catalyst(τ3_hold, k, mixed3_final)
    _, product3_final = split_stream(out3_final, R / (1 + R))

    closed_ext = pfr_step_growth_external_catalyst_recycle(τ3, R, k, c_fresh)
    p_ext_general = 1 - product3_final.x.c / c_fresh
    @test isapprox(p_ext_general, closed_ext.p; rtol=1e-6)

    # --- Step-growth, self-catalyzed ---
    function loop_self(recycle)
        mixed = mix(s_fresh3, recycle)
        out = pfr_step_growth_self_catalyzed(τ3_hold, k, mixed)
        new_recycle, _ = split_stream(out, R / (1 + R))
        return new_recycle
    end
    res_self = solve_tear(loop_self, Stream(R * Q_fresh, s_fresh3.x))
    @test res_self.converged

    mixed4_final = mix(s_fresh3, res_self.stream)
    out4_final = pfr_step_growth_self_catalyzed(τ3_hold, k, mixed4_final)
    _, product4_final = split_stream(out4_final, R / (1 + R))

    closed_self = pfr_step_growth_self_catalyzed_recycle(τ3, R, k, c_fresh)
    p_self_general = 1 - product4_final.x.c / c_fresh
    @test isapprox(p_self_general, closed_self.p; rtol=1e-6)
end

@testset "flowsheet: topology beyond a single recycle loop" begin
    # CSTR -> splitter -> two parallel PFRs at different residence times ->
    # remixed. No recycle here, so this checks composability itself (chaining
    # Stream-based unit operations through a genuine branch-and-remix
    # topology) rather than solve_tear.
    k_p, k_d, k_t, f = 1e3, 1e-5, 1e7, 0.5
    s_feed = Stream(1.0, (I=0.01, M=5.0))

    after_cstr = cstr_free_radical(1000.0, k_p, k_d, k_t, f, s_feed)
    branch1, branch2 = split_stream(after_cstr, 0.5)
    out1 = pfr_free_radical(500.0, k_p, k_d, k_t, f, branch1)
    out2 = pfr_free_radical(2000.0, k_p, k_d, k_t, f, branch2)
    combined = mix(out1, out2)

    # Splitting preserves composition, so both branches start from the same
    # composition; a longer residence time converts more, so the branches'
    # outlet M's must differ, and the remixed M must be their flow-weighted
    # average (mix's own definition) and therefore lie between them.
    @test branch1.x.M == branch2.x.M == after_cstr.x.M
    @test out1.x.M != out2.x.M
    expected_M = (out1.Q * out1.x.M + out2.Q * out2.x.M) / (out1.Q + out2.Q)
    @test isapprox(combined.x.M, expected_M; atol=1e-12)
    @test min(out1.x.M, out2.x.M) <= combined.x.M <= max(out1.x.M, out2.x.M)
    @test isapprox(combined.Q, s_feed.Q; atol=1e-12)  # constant-density mass balance
end

@testset "sanchez-lacombe: constant-hole-volume mixture chemical potentials/activities" begin
    poly = species("polystyrene")
    solv = species("toluene")
    T, P = 298.15, 0.1

    # Pure-component limits: every quantity feeding the chemical potential
    # (Tstar, Pstar, r, v0, x, y) reduces exactly to the pure-component value
    # at w1=1/w1=0 by construction, so ln(a) of the *present* component must
    # be exactly 0 there (the absent component's ln(a) diverges to -Inf,
    # mirroring log(phi2) -> -Inf in ln_activity_polymer as phi2 -> 0).
    r1 = sl_mixture_ln_activities(T, P, solv, 1.0, poly, 0.0)
    @test isapprox(r1.ln_a1, 0.0; atol=1e-9)
    @test r1.ln_a2 == -Inf

    r0 = sl_mixture_ln_activities(T, P, solv, 0.0, poly, 1.0)
    @test isapprox(r0.ln_a2, 0.0; atol=1e-9)
    @test r0.ln_a1 == -Inf

    # Intermediate compositions: finite/real ln_a1 (a2 legitimately
    # underflows toward exp(-thousands) for a long polymer chain -- this is
    # the same qualitative behavior as ln_activity_polymer's "per mole of
    # chain" convention in flory_huggins.jl, not a bug).
    for w1 in (0.05, 0.2, 0.5, 0.8, 0.95)
        r = sl_mixture_ln_activities(T, P, solv, w1, poly, 1 - w1)
        @test isfinite(r.ln_a1)
        @test !isnan(r.ln_a2)
        a = sl_mixture_activities(T, P, solv, w1, poly, 1 - w1)
        @test a.a1 >= 0 && a.a2 >= 0
    end

    # k12 must have a real effect on the activities at fixed composition.
    r_k0 = sl_mixture_ln_activities(T, P, solv, 0.5, poly, 0.5; k12=0.0)
    r_k1 = sl_mixture_ln_activities(T, P, solv, 0.5, poly, 0.5; k12=0.02)
    @test !isapprox(r_k0.ln_a1, r_k1.ln_a1; rtol=1e-6)

    # The central thermodynamic-consistency check this module exists to
    # satisfy (the source's own definition of "consistent"): the Euler
    # relation Sigma_k n_k mu_k - PV = F must hold *exactly* (to floating
    # point) whenever mu_k and P are genuinely partial derivatives of the
    # same free energy F. Checked directly against the actual code (not
    # just derived by hand) at several different states, using the same
    # free-energy expression (Eq. 13 of von Konigslow, Park & Thompson,
    # Phys. Rev. Applied 8, 044009 (2017)) the module's docstring derives
    # mu_k from.
    R_GAS = 8.314462618
    function euler_residual(T, P, sp1, w1, sp2, w2, k12, n1, n2)
        mix = PolyRigorous.sl_consistent_mixing_rules(sp1, w1, sp2, w2; k12=k12)
        ρ̃_mix = PolyRigorous._sl_reduced_density_explicit_r(T, P, mix.Tstar, mix.Pstar, mix.r)

        r1_, r2_ = segment_number(sp1), segment_number(sp2)
        ν1, ν2 = PolyRigorous._sl_segment_volume(sp1), PolyRigorous._sl_segment_volume(sp2)
        x1_charvol, x2_charvol = r1_ * ν1, r2_ * ν2
        Vstar = n1 * x1_charvol + n2 * x2_charvol
        V = Vstar / ρ̃_mix
        n0 = (V - Vstar) / mix.v0

        φ1 = n1 * x1_charvol / V
        φ2 = n2 * x2_charvol / V

        ε11_over_RT = sp1.Tstar / T
        ε22_over_RT = sp2.Tstar / T
        ε12_over_RT = sqrt(sp1.Tstar * sp2.Tstar) * (1 - k12) / T

        F_over_RT = -(V / mix.v0) * (ε11_over_RT * φ1^2 + 2 * ε12_over_RT * φ1 * φ2 + ε22_over_RT * φ2^2) +
                    (n0 * log(n0 * mix.v0 * ℯ / V) - n0) +
                    (n1 * log(φ1) + n2 * log(φ2))

        α1, α2 = r1_ * ν1 / mix.v0, r2_ * ν2 / mix.v0
        μ1_over_RT = PolyRigorous._sl_chem_pot_over_RT(α1, φ1, 1 - ρ̃_mix, ε11_over_RT, ε12_over_RT, φ2)
        μ2_over_RT = PolyRigorous._sl_chem_pot_over_RT(α2, φ2, 1 - ρ̃_mix, ε22_over_RT, ε12_over_RT, φ1)

        P_over_RT = P * 1e6 / (R_GAS * T)
        return (n1 * μ1_over_RT + n2 * μ2_over_RT - P_over_RT * V) - F_over_RT
    end

    for (Tc, Pc, w1c, k12c) in [(298.15, 0.1, 0.4, 0.0), (350.0, 5.0, 0.15, 0.02), (280.0, 20.0, 0.7, -0.01)]
        n1c, n2c = 1000 * (w1c / solv.M), 1000 * ((1 - w1c) / poly.M)
        resid = euler_residual(Tc, Pc, solv, w1c, poly, 1 - w1c, k12c, n1c, n2c)
        @test isapprox(resid, 0.0; atol=1e-6)
    end
end

@testset "reactors: staged initiator injection (LDPE autoclave/tubular)" begin
    k_p, k_d, k_t, f = 1e3, 1e-5, 1e7, 0.5
    I0, M0 = 0.01, 5.0
    τ = 3600.0

    # A single stage with the whole charge injected up front must reduce
    # exactly to the plain (non-staged) CSTR/PFR.
    staged1_cstr = cstr_train_free_radical_staged([τ], k_p, k_d, k_t, f, [I0], M0)
    plain_cstr = cstr_free_radical(τ, k_p, k_d, k_t, f, I0, M0)
    @test isapprox(staged1_cstr.M, plain_cstr.M; rtol=1e-12)
    @test isapprox(staged1_cstr.I, plain_cstr.I; rtol=1e-12)

    staged1_pfr = pfr_train_free_radical_staged([τ], k_p, k_d, k_t, f, [I0], M0)
    plain_pfr = pfr_free_radical(τ, k_p, k_d, k_t, f, I0, M0)
    @test isapprox(staged1_pfr.M, plain_pfr.M; rtol=1e-12)

    # Splitting a plain PFR into extra segments with zero additional
    # injection must not change the result at all -- this isolates the
    # staging/carry-over machinery from the injection feature itself.
    n = 7
    τs = fill(τ / n, n)
    ΔI_zero_after_first = vcat([I0], zeros(n - 1))
    staged_n = pfr_train_free_radical_staged(τs, k_p, k_d, k_t, f, ΔI_zero_after_first, M0)
    @test isapprox(staged_n.M, plain_pfr.M; rtol=1e-8)
    @test isapprox(staged_n.I, I0 * exp(-k_d * τ); rtol=1e-8)

    # Injecting strictly more total initiator (an extra pulse on top of the
    # same first-stage charge) must not decrease conversion, for both the
    # autoclave (CSTR-train) and tubular (PFR-train) models.
    n2 = 4
    τs2 = fill(τ / n2, n2)
    ΔI_baseline = vcat([I0], zeros(n2 - 1))
    ΔI_extra = vcat([I0], [I0], zeros(n2 - 2))
    @test cstr_train_free_radical_staged(τs2, k_p, k_d, k_t, f, ΔI_extra, M0).conversion >
          cstr_train_free_radical_staged(τs2, k_p, k_d, k_t, f, ΔI_baseline, M0).conversion
    @test pfr_train_free_radical_staged(τs2, k_p, k_d, k_t, f, ΔI_extra, M0).conversion >
          pfr_train_free_radical_staged(τs2, k_p, k_d, k_t, f, ΔI_baseline, M0).conversion

    @test_throws ArgumentError cstr_train_free_radical_staged([τ, τ], k_p, k_d, k_t, f, [I0], M0)
    @test_throws ArgumentError pfr_train_free_radical_staged([τ, τ], k_p, k_d, k_t, f, [I0], M0)

    # Ideal-gas concentration helper: [M] = P/(RT), a basic dimensional
    # check against the gas constant directly (1 mol of ideal gas at its
    # own characteristic P,T satisfies PV=RT, i.e. concentration 1/V=P/(RT)).
    R_GAS = 8.314462618
    T = 350.0
    P = 2.0  # MPa
    c = ideal_gas_concentration(P, T)
    @test isapprox(c, P * 1e6 / (R_GAS * T) / 1000; rtol=1e-12)
    @test c > 0
    # Doubling pressure at fixed T must exactly double the concentration.
    @test isapprox(ideal_gas_concentration(2P, T), 2c; rtol=1e-12)
end

end # @testset "PolyRigorous"
