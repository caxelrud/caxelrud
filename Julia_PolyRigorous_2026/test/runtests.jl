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

end # @testset "PolyRigorous"
