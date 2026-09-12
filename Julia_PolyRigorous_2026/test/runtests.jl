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
    # (Raoult's-law-like) behavior ln(a1) = ln(phi1).
    phi2 = 0.3
    @test isapprox(PolyRigorous.ln_activity_solvent(phi2, 1.0, 0.0), log(1 - phi2); atol=1e-12)

    # As N -> infinity at fixed phi2, chi (the (1-1/N) term saturates)
    ln_a1_bigN = PolyRigorous.ln_activity_solvent(phi2, 1e8, 0.5)
    ln_a1_inftyN = log(1 - phi2) + phi2 + 0.5 * phi2^2
    @test isapprox(ln_a1_bigN, ln_a1_inftyN; atol=1e-6)

    # Symmetry check for N=1: solvent and polymer are the same size, so
    # relabeling which species is "solvent" (phi2 <-> 1-phi2) should map
    # ln_activity_solvent onto ln_activity_polymer.
    @test isapprox(
        PolyRigorous.ln_activity_solvent(phi2, 1.0, 0.7),
        PolyRigorous.ln_activity_polymer(1 - phi2, 1.0, 0.7);
        atol=1e-12,
    )
end

@testset "flory-huggins: critical point and spinodal" begin
    N = 1.0
    crit = critical_point(N)
    @test isapprox(crit.phi2c, 0.5; atol=1e-12)
    @test isapprox(crit.chic, 2.0; atol=1e-12)

    N2 = 1000.0
    crit2 = critical_point(N2)
    @test isapprox(crit2.phi2c, 1 / (1 + sqrt(N2)); atol=1e-12)

    # The spinodal touches the binodal exactly at the critical point:
    # chi_spinodal(phi2c) == chic.
    @test isapprox(chi_spinodal(crit2.phi2c, N2), crit2.chic; atol=1e-9)
end

@testset "phase equilibrium: binodal brackets the critical point" begin
    N = 50.0
    crit = critical_point(N)
    chi = 1.05 * crit.chic
    res = binodal_pair(N, chi)
    @test res.converged
    @test res.phi2a < crit.phi2c < res.phi2b

    # Both branches should satisfy the equal-activity conditions to high
    # precision.
    a1a = PolyRigorous.ln_activity_solvent(res.phi2a, N, chi)
    a1b = PolyRigorous.ln_activity_solvent(res.phi2b, N, chi)
    @test isapprox(a1a, a1b; atol=1e-6)

    a2a = PolyRigorous.ln_activity_polymer(res.phi2a, N, chi)
    a2b = PolyRigorous.ln_activity_polymer(res.phi2b, N, chi)
    @test isapprox(a2a, a2b; atol=1e-6)
end

@testset "phase equilibrium: binodal_curve and spinodal_curve run and stay ordered" begin
    N = 80.0
    curve = binodal_curve(N; npoints=15)
    @test length(curve) > 5
    for pt in curve
        @test pt.phi2a < pt.phi2b
    end

    spin = spinodal_curve(N; npoints=50)
    @test length(spin) == 50
    @test all(pt.chi > 0 for pt in spin)
end

@testset "sanchez-lacombe: EOS residual and PVT sanity" begin
    tol = species("toluene")
    T, P = 298.15, 0.1  # K, MPa (~1 atm)

    rho_t = PolyRigorous.reduced_density(T, P, tol)
    @test 0 < rho_t < 1

    r = segment_number(tol)
    resid = sl_eos_residual(rho_t, T / tol.Tstar, P / tol.Pstar, r)
    @test isapprox(resid, 0.0; atol=1e-8)

    rho = density(T, P, tol)
    @test 0 < rho < tol.rhostar

    # Increasing pressure at fixed T should compress the fluid (density
    # goes up).
    rho_hi = density(T, 50.0, tol)
    @test rho_hi > rho

    # Increasing temperature at fixed P should expand the fluid (density
    # goes down) below the polymer's/solvent's spinodal-like limits.
    rho_hot = density(T + 20, P, tol)
    @test rho_hot < rho

    alpha = thermal_expansion_coefficient(T, P, tol)
    @test alpha > 0

    kappa = isothermal_compressibility(T, P, tol)
    @test kappa > 0
end

@testset "free-radical kinetics" begin
    kp, kd, kt, f = 1e3, 1e-5, 1e7, 0.5
    I0, M0 = 0.01, 5.0

    # No reaction has happened yet at t=0.
    @test monomer_concentration(0.0, kp, kd, kt, f, I0, M0) == M0
    @test conversion(0.0, kp, kd, kt, f, I0, M0) == 0.0

    # Conversion increases monotonically in t and never exceeds 1 -- and, in
    # this dead-end (finite initiator charge) batch model, never reaches it
    # either: once the initiator is exhausted, propagation stops with some
    # monomer left unreacted. The limiting conversion as t -> infinity has a
    # known closed form (exp(-kd t/2) -> 0 in the conversion formula).
    ts = [0.0, 100.0, 1e4, 1e6, 1e9]
    convs = [conversion(t, kp, kd, kt, f, I0, M0) for t in ts]
    @test issorted(convs)
    @test all(0 .<= convs .< 1)
    conv_limit = 1 - exp(-(2 * kp / kd) * sqrt(f * kd * I0 / kt))
    @test isapprox(convs[end], conv_limit; atol=1e-9)

    # Two independent routes to the kinetic chain length should agree:
    # nu = Rp/Ri directly, vs. the closed-form kinetic_chain_length formula.
    Mrad = radical_concentration(f, kd, I0, kt)
    Ri = initiation_rate(f, kd, I0)
    Rp = propagation_rate(kp, M0, Mrad)
    nu_direct = Rp / Ri
    nu_formula = kinetic_chain_length(kp, M0, f, kd, I0, kt)
    @test isapprox(nu_direct, nu_formula; rtol=1e-10)

    # Xn_mixed reduces to the pure combination/disproportionation cases at
    # its endpoints.
    nu = nu_formula
    @test isapprox(Xn_mixed(nu, 0.0), Xn_combination(nu); atol=1e-10)
    @test isapprox(Xn_mixed(nu, 1.0), Xn_disproportionation(nu); atol=1e-10)
    @test Xn_combination(nu) > Xn_disproportionation(nu)  # 2nu > nu
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

end # @testset "PolyRigorous"
