using Test
using PolyRigorous

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

    # Symmetry check for N=1: swapping the roles of "solvent" and "polymer"
    # should give the same activity expression.
    @test isapprox(
        PolyRigorous.ln_activity_solvent(phi2, 1.0, 0.7),
        PolyRigorous.ln_activity_polymer(phi2, 1.0, 0.7);
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
