### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 1d69aa37-9d2a-4af7-a1d3-9c791a983c12
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    using PlutoUI
    using Plots
    using PolyRigorous
end

# ╔═╡ ba17fea8-391c-4a6c-9657-ca1a468b3047
md"""
# PolyRigorous — Thermodynamics Explorer

Interactive exploration of the Flory-Huggins and Sanchez-Lacombe polymer
solution thermodynamics implemented in the `PolyRigorous` package
(`../src`). This notebook is a thin UI on top of that package — all the
physics lives in the package, not here, so it stays testable independent
of Pluto.
"""

# ╔═╡ bf44a61f-3c6f-4212-8320-922d3dd222dc
md"## 1. Choose a polymer / solvent pair"

# ╔═╡ 8546f5e2-21e0-4173-bbad-67fe6a4b4c8d
begin
    polymer_options = [k => SPECIES_DB[k].name for k in sort(collect(keys(SPECIES_DB))) if SPECIES_DB[k].kind == :polymer]
    solvent_options = [k => SPECIES_DB[k].name for k in sort(collect(keys(SPECIES_DB))) if SPECIES_DB[k].kind == :solvent]
end

# ╔═╡ d8a30bc1-4969-4055-9d55-1efc5369b5d6
@bind polymer_key Select(polymer_options; default="polystyrene")

# ╔═╡ e4e85e1c-02b7-4272-ab76-ba988ccda7d3
@bind solvent_key Select(solvent_options; default="toluene")

# ╔═╡ fd4a18c4-1ad2-466e-b889-98b2c8a5ff68
@bind T_C Slider(0:1:150; default=25, show_value=true)

# ╔═╡ 7651a9e8-d0c8-435f-8e1d-21dc48be24dd
begin
    poly = species(polymer_key)
    solv = species(solvent_key)
    N = degree_of_polymerization(poly, solv)
    T = T_C + 273.15
end

# ╔═╡ d9924eaa-de0e-43ee-b9af-eeba20ea4ff2
md"""
Selected pair: **$(poly.name)** in **$(solv.name)**, at $(T_C) °C ($(round(T, digits=1)) K).

Flory-Huggins size ratio N = Vm(polymer)/Vm(solvent) = **$(round(N, digits=1))**
"""

# ╔═╡ 68459002-17d6-4071-89f3-997d3f453ce8
md"## 2. Flory-Huggins interaction parameter χ"

# ╔═╡ 664d29e8-9bac-414d-91fb-637e798beb72
@bind χₛ Slider(0.20:0.01:0.50; default=0.34, show_value=true)

# ╔═╡ 068f474e-d08c-4204-b6e9-dd61fd2147c0
begin
    χ = chi_from_solubility(solv.Vm, solv.δ, poly.δ, T; χₛ=χₛ)
    crit = critical_point(N)
end

# ╔═╡ 1c1f35f0-9954-40f5-86aa-05fc5da71c63
md"""
χ (from solubility parameters, entropic correction χₛ = $(χₛ)) = **$(round(χ, digits=3))**

Critical point for this N: φ₂,c = $(round(crit.φ₂_c, digits=4)), χ_c = $(round(crit.χ_c, digits=3))

$(χ > crit.χ_c ? "⚠️ χ > χ_c at this temperature: the model predicts liquid-liquid phase separation over part of the composition range." : "✅ χ < χ_c at this temperature: the model predicts complete miscibility at all compositions.")
"""

# ╔═╡ 715647ac-d886-47bd-ac78-710f2b0ffc41
begin
    φ₂s = range(0.001, 0.999; length=300)
    ln_a1 = [ln_activity_solvent(p, N, χ) for p in φ₂s]
    plot(φ₂s, ln_a1;
        xlabel="polymer volume fraction φ₂", ylabel="ln(a₁)",
        label="ln(solvent activity)", lw=2, legend=:bottomleft,
        title="Solvent activity vs. composition")
end

# ╔═╡ 338160d0-296d-46f3-b760-8f4b2190af0f
md"## 3. Phase diagram (χ vs. φ₂)"

# ╔═╡ ab4ee0c7-8a80-4431-81b4-30ef03130da6
begin
    spin = spinodal_curve(N; npoints=300)
    plt = plot([s.φ₂ for s in spin], [s.χ for s in spin];
        label="spinodal", lw=2, xlabel="φ₂", ylabel="χ",
        title="Phase diagram at N = $(round(N, digits=1))", legend=:topright)
    scatter!(plt, [crit.φ₂_c], [crit.χ_c]; label="critical point", ms=6)
    if χ > crit.χ_c
        bcurve = binodal_curve(N; χ_max=max(1.001χ, 1.2crit.χ_c), npoints=40)
        if !isempty(bcurve)
            plot!(plt, [b.φ₂a for b in bcurve], [b.χ for b in bcurve];
                label="binodal", lw=2, ls=:dash, color=3)
            plot!(plt, [b.φ₂b for b in bcurve], [b.χ for b in bcurve];
                label=nothing, lw=2, ls=:dash, color=3)
        end
    end
    hline!(plt, [χ]; label="current χ", ls=:dot, color=:black)
    plt
end

# ╔═╡ 6103bceb-8a74-49ff-8e03-54a0220ad158
md"## 4. Sanchez-Lacombe pure-component PVT"

# ╔═╡ 66b69bd7-af4e-4465-a0dc-6f167b696b18
sl_options = [k => SPECIES_DB[k].name for k in sort(collect(keys(SPECIES_DB)))]

# ╔═╡ 0600f577-7093-4143-8b69-321c0b2fbd31
@bind sl_key Select(sl_options; default="toluene")

# ╔═╡ 045bea12-4120-4a56-afc7-b386a98d8936
@bind P_MPa Slider(0.1:0.5:50; default=0.1, show_value=true)

# ╔═╡ 10f6793d-7f12-4f7a-8b39-67a8d4914c7f
begin
    sl_target = species(sl_key)
    Ts = range(0.55sl_target.Tstar, 0.95sl_target.Tstar; length=60)
    rhos = Float64[]
    Ts_ok = Float64[]
    for Tv in Ts
        try
            push!(rhos, density(Tv, P_MPa, sl_target))
            push!(Ts_ok, Tv)
        catch
        end
    end
    plot(Ts_ok, rhos;
        xlabel="T (K)", ylabel="density (kg/m³)",
        label="$(sl_target.name) at P = $(P_MPa) MPa", lw=2,
        title="Sanchez-Lacombe density vs. temperature")
end

# ╔═╡ 4b222d62-b1d1-4c02-9f3a-85461ce59b1c
md"""
## 5. Free-radical polymerization kinetics

Batch, isothermal free-radical polymerization under the standard
quasi-steady-state approximation (QSSA), from `kinetics_free_radical.jl`.
Propagation (`k_p`) and termination (`k_t`) rate constants are held fixed
at typical styrene-like values; vary the initiator decomposition rate
`k_d`, initial initiator concentration `[I]₀`, and the disproportionation
fraction `δ` of termination events.
"""

# ╔═╡ fe3637be-868c-4ee2-bc79-899656be2deb
@bind log10_k_d Slider(-6:0.25:-3; default=-5, show_value=true)

# ╔═╡ 151fb00f-1d16-47fc-b605-8f5e095d1665
@bind I0_frk Slider(0.001:0.001:0.05; default=0.01, show_value=true)

# ╔═╡ 0ad80cab-531b-405d-8837-b9074d32fcc4
@bind δ_term Slider(0.0:0.05:1.0; default=0.2, show_value=true)

# ╔═╡ 2689f077-50ab-4917-9797-2b3c8bce3551
begin
    k_p_frk = 1.0e3  # L/(mol s), typical propagation rate constant
    k_t_frk = 1.0e7  # L/(mol s), typical termination rate constant
    f_frk = 0.5      # initiator efficiency
    M0_frk = 5.0     # mol/L, bulk-ish monomer concentration
    k_d_frk = 10.0^log10_k_d
    ν_frk = kinetic_chain_length(k_p_frk, M0_frk, f_frk, k_d_frk, I0_frk, k_t_frk)
    Xn_frk = Xn_mixed(ν_frk, δ_term)
end

# ╔═╡ 7534ac44-9487-48af-883d-6e8805b60d70
md"""
k_d = $(round(k_d_frk, sigdigits=3)) 1/s, [I]₀ = $(I0_frk) mol/L, δ = $(δ_term)

Kinetic chain length ν = $(round(ν_frk, digits=1)); number-average degree of polymerization Xₙ = $(round(Xn_frk, digits=1))

(Instantaneous PDI is exactly 1.5 for pure combination (δ=0) and 2.0 for pure disproportionation (δ=1); this package does not yet provide a closed-form PDI for the intermediate mixed case.)
"""

# ╔═╡ 5d467096-e27f-4b7e-9d9b-4c33d26e590e
begin
    ts_frk = range(0, 5 / k_d_frk; length=300)  # a few initiator half-lives
    convs_frk = [conversion(t, k_p_frk, k_d_frk, k_t_frk, f_frk, I0_frk, M0_frk) for t in ts_frk]
    plot(ts_frk, convs_frk;
        xlabel="time (s)", ylabel="monomer conversion",
        label=nothing, lw=2,
        title="Batch free-radical conversion vs. time")
end

# ╔═╡ 846ed143-a1ae-4fff-b383-e5404b9112bb
md"""
## 6. Step-growth polymerization & the Flory distribution

Extent of reaction `p(t)` for the two classic step-growth kinetic cases
(externally catalyzed = second order overall; self-catalyzed = third
order overall), from `kinetics_step_growth.jl`, plus the resulting Flory
"most probable" molecular weight distribution at the extent of reaction
reached by time `t_max`.
"""

# ╔═╡ e5d1ac62-d0e6-4efa-aa5a-302ea607b969
@bind k_step Slider(0.01:0.01:2.0; default=0.5, show_value=true)

# ╔═╡ 12428167-5594-4efa-bef9-175a32d67f1d
@bind c0_step Slider(0.5:0.5:5.0; default=1.0, show_value=true)

# ╔═╡ 69263b58-861a-41a0-81bd-6718c32ad6fb
@bind t_max_step Slider(1:1:200; default=50, show_value=true)

# ╔═╡ 52849a69-c1ac-418e-8f52-35958b930477
begin
    p_ext_step = extent_reaction_external_catalyst(k_step, c0_step, t_max_step)
    p_self_step = extent_reaction_self_catalyzed(k_step, c0_step, t_max_step)
    Xn_ext_step = carothers_Xn(p_ext_step)
    Xn_self_step = carothers_Xn(p_self_step)
end

# ╔═╡ f2306b92-42be-4d25-a276-55e45bb8401d
md"""
At t = $(t_max_step): p (external catalyst) = $(round(p_ext_step, digits=4)), Xₙ = $(round(Xn_ext_step, digits=1)); p (self-catalyzed) = $(round(p_self_step, digits=4)), Xₙ = $(round(Xn_self_step, digits=1)).

PDI = 1+p → $(round(flory_PDI(p_ext_step), digits=3)) (external catalyst case), approaching 2 as p → 1.
"""

# ╔═╡ 38670f3b-dd5a-4a13-91dd-49db4586e034
begin
    ts_step = range(0, t_max_step; length=300)
    p_ext_curve = [extent_reaction_external_catalyst(k_step, c0_step, t) for t in ts_step]
    p_self_curve = [extent_reaction_self_catalyzed(k_step, c0_step, t) for t in ts_step]
    plot(ts_step, p_ext_curve; label="external catalyst (2nd order)", lw=2,
        xlabel="time", ylabel="extent of reaction p", legend=:bottomright,
        title="Step-growth extent of reaction vs. time")
    plot!(ts_step, p_self_curve; label="self-catalyzed (3rd order)", lw=2, ls=:dash)
end

# ╔═╡ 0028e23f-fa21-4be0-a9f9-51a727765eff
begin
    xmax_flory = min(ceil(Int, 10 * Xn_ext_step), 500)
    xs_flory = 1:xmax_flory
    mole_fracs = [flory_mole_fraction(x, p_ext_step) for x in xs_flory]
    weight_fracs = [flory_weight_fraction(x, p_ext_step) for x in xs_flory]
    plot(xs_flory, mole_fracs; label="mole fraction", lw=2,
        xlabel="chain length x", ylabel="fraction",
        title="Flory most-probable distribution (external-catalyst branch)")
    plot!(xs_flory, weight_fracs; label="weight fraction", lw=2, ls=:dash)
end

# ╔═╡ 065d37c5-ab38-430f-a26b-3682ea251778
md"""
## 7. Reactor models: batch/PFR vs. CSTR

An ideal PFR is kinetically equivalent to a batch reactor run for a time
equal to its residence time τ — both curves above (free-radical
conversion vs. time, step-growth extent of reaction vs. time) are exactly
that PFR/batch curve. A CSTR, fed continuously and mixed instantaneously,
instead runs its *entire* volume at one steady-state (outlet) composition
— so for these kinetics, where the rate falls as the reactant is consumed,
a single CSTR is always less efficient than a PFR/batch at the same
residence time: it spends the whole time reacting at the low outlet
concentration instead of starting fast at the feed concentration. Both
panels below reuse the sliders from sections 5 and 6 — see
`reactors.jl`.
"""

# ╔═╡ f6c82001-9a3d-49ab-871b-3a8c77f26a95
begin
    cstr_convs_frk = [cstr_free_radical(τ, k_p_frk, k_d_frk, k_t_frk, f_frk, I0_frk, M0_frk).conversion for τ in ts_frk]
    plot(ts_frk, convs_frk;
        xlabel="residence time τ (s)", ylabel="conversion", legend=:bottomright,
        label="PFR / batch", lw=2, title="Free-radical: CSTR vs. PFR/batch conversion")
    plot!(ts_frk, cstr_convs_frk; label="CSTR", lw=2, ls=:dash)
end

# ╔═╡ 3ae19407-0464-49c4-a3f1-0b305be274d2
begin
    cstr_p_ext_curve = [cstr_step_growth_external_catalyst(t, k_step, c0_step).p for t in ts_step]
    plot(ts_step, p_ext_curve;
        xlabel="residence time τ", ylabel="extent of reaction p", legend=:bottomright,
        label="PFR / batch (external catalyst)", lw=2,
        title="Step-growth: CSTR vs. PFR/batch extent of reaction")
    plot!(ts_step, cstr_p_ext_curve; label="CSTR (external catalyst)", lw=2, ls=:dash)
end

# ╔═╡ d9972239-b3e5-4eba-bd55-d4f90c770282
md"""
## 8. Coordination polymerization kinetics

Ziegler-Natta / metallocene-type coordination polymerization, from
`kinetics_coordination.jl`. Unlike free-radical termination, chain
transfer here doesn't kill the active site — it releases one dead chain
and immediately starts a new one — so Xₙ is set by the ratio of
propagation to the *total* rate of chain-releasing events (transfer +
true termination), not a bimolecular termination step. Vary the active
site concentration [C*]₀, the true-termination rate k_t (0 = fully
"living", no site deactivation), and the hydrogen concentration [H₂] —
industrially, H₂ is dosed specifically to *lower* molecular weight via
chain transfer to hydrogen.
"""

# ╔═╡ 7ad99fa2-83f8-450f-8f70-e4dcc2c6edfb
@bind log10_C0_star Slider(-6:0.25:-3; default=-4, show_value=true)

# ╔═╡ 48b3199c-d5c5-44a6-8ef6-fd93e59d4e78
@bind k_t_coord Slider(0.0:0.0002:0.01; default=0.001, show_value=true)

# ╔═╡ 6f8f3320-8e35-4c43-b7c4-6f5d2d4ff11f
@bind H2_coord Slider(0.0:0.01:1.0; default=0.1, show_value=true)

# ╔═╡ 1e0ff82a-d49e-4c10-bf22-a0fc4edcd9fb
begin
    k_p_coord = 50.0      # L/(mol s), illustrative propagation rate constant
    k_trM_coord = 0.05    # L/(mol s), transfer to monomer
    k_trH_coord = 0.3     # L/(mol s), transfer to hydrogen
    k_tr0_coord = 1.0e-3  # 1/s, spontaneous transfer
    M0_coord = 5.0        # mol/L, bulk-ish monomer concentration
    C0_star_coord = 10.0^log10_C0_star
    k_release_coord = transfer_rate_constant(;
        k_trM=k_trM_coord, M=M0_coord, k_trH=k_trH_coord, H2=H2_coord,
        k_tr0=k_tr0_coord, k_t=k_t_coord,
    )
    Xn_coord = Xn_coordination(k_p_coord, M0_coord, k_release_coord)
end

# ╔═╡ 56ba9de1-48e7-4252-9b7d-dc4a1e7900dd
md"""
[C*]₀ = $(round(C0_star_coord, sigdigits=3)) mol/L, k_t = $(round(k_t_coord, sigdigits=3)) 1/s, [H₂] = $(H2_coord) mol/L

Instantaneous Xₙ = $(round(Xn_coord, digits=1)) (at the initial monomer concentration)

$(k_t_coord == 0 ? "✅ k_t = 0: a fully \"living\" system — the active-site pool never decays, so given enough time all the monomer is eventually consumed." : "The active-site pool decays with a half-life of $(round(log(2)/k_t_coord, sigdigits=3)) s, so conversion will plateau below 100% (a \"dead-end\" polymerization, like the free-radical case in section 5).")
"""

# ╔═╡ da08bb91-9401-4b31-9952-ba1fccf997ff
begin
    ts_coord = k_t_coord == 0 ? range(0, 5 / (k_p_coord * C0_star_coord); length=300) : range(0, 5 / k_t_coord; length=300)
    convs_coord = [conversion_coordination(t, k_p_coord, C0_star_coord, k_t_coord, M0_coord) for t in ts_coord]
    plot(ts_coord, convs_coord;
        xlabel="time (s)", ylabel="monomer conversion", label=nothing, lw=2,
        title="Batch coordination polymerization conversion vs. time")
end

# ╔═╡ 4bd841d5-5e7d-48c8-994a-2f48bc049dc5
begin
    H2_range = range(0.0, 2.0; length=200)
    Xn_vs_H2 = [Xn_coordination(k_p_coord, M0_coord, transfer_rate_constant(;
        k_trM=k_trM_coord, M=M0_coord, k_trH=k_trH_coord, H2=h,
        k_tr0=k_tr0_coord, k_t=k_t_coord,
    )) for h in H2_range]
    plot(H2_range, Xn_vs_H2;
        xlabel="[H₂] (mol/L)", ylabel="instantaneous Xₙ", label=nothing, lw=2,
        title="Hydrogen response: Xₙ vs. [H₂] (industrial MW control lever)")
end

# ╔═╡ 8139da23-ae9c-4052-8bd3-32ef4fa0391f
md"""
## 9. Reactor trains: CSTRs in series → PFR

Chaining `n` ideal CSTRs in series, each stage's outlet feeding the next
stage's inlet (`reactors.jl`), and holding the *total* residence time
fixed while increasing `n`: performance should climb monotonically from
the single-CSTR value toward the PFR/batch value as `n → ∞` — a classic
reactor-engineering result. Reuses the free-radical parameters from
section 5, with total residence time equal to that section's time axis.
"""

# ╔═╡ df12ac0e-1ac5-4ed5-80bc-78be20d8940c
@bind n_stages_frk Slider(1:1:50; default=5, show_value=true)

# ╔═╡ fefabb16-4409-48a6-b6d6-10ae48e15e97
begin
    τ_total_frk = ts_frk[end]
    train_res_frk = cstr_train_free_radical(n_stages_frk, τ_total_frk / n_stages_frk, k_p_frk, k_d_frk, k_t_frk, f_frk, I0_frk, M0_frk)
    single_cstr_frk = cstr_free_radical(τ_total_frk, k_p_frk, k_d_frk, k_t_frk, f_frk, I0_frk, M0_frk)
    pfr_frk = pfr_free_radical(τ_total_frk, k_p_frk, k_d_frk, k_t_frk, f_frk, I0_frk, M0_frk)
end

# ╔═╡ 142995ab-27e7-45fe-9cd4-2e7a6bd60151
md"""
At fixed total residence time τ = $(round(τ_total_frk, sigdigits=3)) s: a single CSTR gives conversion = $(round(single_cstr_frk.conversion, digits=4)); $(n_stages_frk) equal CSTRs in series give $(round(train_res_frk.conversion, digits=4)); the PFR/batch limit is $(round(pfr_frk.conversion, digits=4)).
"""

# ╔═╡ 91358a8b-c845-4c13-9779-d9d01ac95cd6
begin
    ns_frk = 1:50
    train_convs_frk = [cstr_train_free_radical(n, τ_total_frk / n, k_p_frk, k_d_frk, k_t_frk, f_frk, I0_frk, M0_frk).conversion for n in ns_frk]
    plot(ns_frk, train_convs_frk;
        xlabel="number of equal CSTR stages", ylabel="conversion", label="CSTR train", lw=2,
        legend=:bottomright, title="CSTR train → PFR as stage count grows")
    hline!([pfr_frk.conversion]; label="PFR / batch limit", ls=:dash, color=:black)
end

# ╔═╡ 0111a035-58a0-45ed-8913-0fea5774366e
md"""
## 10. Sanchez-Lacombe binary mixture PVT

Mixture density from the Sanchez-Lacombe mixing rules (`sl_mixing_rules`,
`sl_mixture_species`, `sanchez_lacombe_mixture.jl`), applied to the same
polymer/solvent pair and temperature chosen in section 1 and the pressure
from section 4. `w₁` is the *weight fraction of solvent* in the mixture.
Only PVT (density) behavior is covered here; mixture chemical
potentials/activities use a *different* (thermodynamically consistent)
mixing-rule convention — see section 13.
"""

# ╔═╡ bba86433-0285-4076-8b58-58e6c870362d
@bind w1_mix Slider(0.0:0.02:1.0; default=0.5, show_value=true)

# ╔═╡ 3977b7bf-5ee5-49b0-a831-9a9332e7e73d
begin
    sp_mix = sl_mixture_species(solv, w1_mix, poly, 1 - w1_mix)
    ρ_mix = density(T, P_MPa, sp_mix)
end

# ╔═╡ 9497e799-3098-446c-892e-b485ac224e9a
md"""
$(round(100w1_mix, digits=0))% $(solv.name) / $(round(100 * (1 - w1_mix), digits=0))% $(poly.name) by mass, at $(T_C) °C, P = $(P_MPa) MPa: mixture density = $(round(ρ_mix, digits=1)) kg/m³ (pure $(solv.name): $(round(density(T, P_MPa, solv), digits=1)) kg/m³; pure $(poly.name): $(round(density(T, P_MPa, poly), digits=1)) kg/m³).
"""

# ╔═╡ 6712044c-4a44-4e58-b3a9-d00604247b2f
begin
    w1s = range(0.0, 1.0; length=100)
    ρs_mix = [density(T, P_MPa, sl_mixture_species(solv, w, poly, 1 - w)) for w in w1s]
    plot(w1s, ρs_mix;
        xlabel="w₁ (solvent weight fraction)", ylabel="density (kg/m³)",
        label=nothing, lw=2,
        title="Sanchez-Lacombe mixture density vs. composition")
end

# ╔═╡ 73fbbb0f-2400-43e9-8e9e-4a9161a21414
md"""
## 11. PFR with recycle

A single ideal PFR whose outlet stream is split, with a fraction
`R = Q_recycle/Q_fresh` sent back and mixed into the fresh feed
(`reactor_recycle.jl`) — the first step beyond a single pass toward
genuine flowsheet topology, since the recycle stream's composition depends
on the very reactor outlet it feeds, making this an implicit
("flowsheet-style") problem. At `R = 0` this is exactly the plain PFR
curve from sections 7/9 above; as `R → ∞` it approaches the CSTR at the
*same* nominal τ — a classic reactor-engineering result (a PFR with
infinite recycle behaves like a CSTR), checked directly in the test suite.
Reuses the free-radical parameters and total residence time from sections
5 and 9, the coordination parameters from section 8, and the step-growth
external-catalyst parameters from section 6.
"""

# ╔═╡ 8822e6f6-0f04-420c-9387-70607c8e6cf4
@bind R_recycle Slider([0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 500.0, 1000.0, 10000.0]; default=1.0, show_value=true)

# ╔═╡ 569b1802-d2fc-4092-a694-11b6703fc640
begin
    τ_recycle_coord = 100.0  # s, illustrative fixed residence time for the coordination case
    τ_recycle_step = t_max_step

    rec_frk = pfr_free_radical_recycle(τ_total_frk, R_recycle, k_p_frk, k_d_frk, k_t_frk, f_frk, I0_frk, M0_frk)
    rec_coord = pfr_coordination_recycle(τ_recycle_coord, R_recycle, k_p_coord, k_t_coord, C0_star_coord, M0_coord)
    rec_step = pfr_step_growth_external_catalyst_recycle(τ_recycle_step, R_recycle, k_step, c0_step)

    pfr_limit_coord = pfr_coordination(τ_recycle_coord, k_p_coord, C0_star_coord, k_t_coord, M0_coord)
    cstr_limit_coord = cstr_coordination(τ_recycle_coord, k_p_coord, k_t_coord, C0_star_coord, M0_coord)
    pfr_limit_step = pfr_step_growth_external_catalyst(τ_recycle_step, k_step, c0_step)
    cstr_limit_step = cstr_step_growth_external_catalyst(τ_recycle_step, k_step, c0_step)
end

# ╔═╡ dd2436b1-9501-4982-9d6b-3f3649c45100
md"""
At R = $(R_recycle):

- Free-radical (τ = $(round(τ_total_frk, sigdigits=3)) s): conversion = $(round(rec_frk.conversion, digits=4)) (PFR limit at R=0: $(round(pfr_frk.conversion, digits=4)); CSTR limit at R→∞: $(round(single_cstr_frk.conversion, digits=4)))
- Coordination (τ = $(τ_recycle_coord) s): conversion = $(round(rec_coord.conversion, digits=4)) (PFR limit: $(round(pfr_limit_coord.conversion, digits=4)); CSTR limit: $(round(cstr_limit_coord.conversion, digits=4)))
- Step-growth, external catalyst (τ = $(τ_recycle_step)): p = $(round(rec_step.p, digits=4)) (PFR limit: $(round(pfr_limit_step.p, digits=4)); CSTR limit: $(round(cstr_limit_step.p, digits=4)))
"""

# ╔═╡ 66ee83fc-148a-4a85-91a0-fbba76f70197
begin
    Rs_plot = exp10.(range(-2, 3; length=100))
    convs_recycle_frk = [pfr_free_radical_recycle(τ_total_frk, R, k_p_frk, k_d_frk, k_t_frk, f_frk, I0_frk, M0_frk).conversion for R in Rs_plot]
    plot(Rs_plot, convs_recycle_frk;
        xscale=:log10, xlabel="recycle ratio R", ylabel="conversion", legend=:right,
        label="PFR with recycle", lw=2, title="Free-radical: PFR-with-recycle conversion vs. R")
    hline!([pfr_frk.conversion]; label="PFR (R=0) limit", ls=:dash, color=:black)
    hline!([single_cstr_frk.conversion]; label="CSTR (R→∞) limit", ls=:dot, color=:red)
end

# ╔═╡ 178d60eb-0691-4815-8822-d1d98ce1a893
md"""
## 12. General flowsheet: branch-and-remix

Beyond a single recycle loop (section 11), `flowsheet.jl`'s `Stream`,
`mix`, and `split_stream` primitives compose into *any* topology — this
one doesn't even need `solve_tear`, since there's no cycle to converge.
The free-radical fresh feed (section 5) is split into two branches, each
run through its own PFR at a different residence time, then remixed —
plain Julia function composition over `Stream`s, reusing the same
[`pfr_free_radical`](@ref) reactor function as everywhere else (via its
`Stream`-taking method).
"""

# ╔═╡ a2a2b30c-5670-48f4-9ba9-54fcd576e200
@bind frac_branch1 Slider(0.0:0.02:1.0; default=0.5, show_value=true)

# ╔═╡ 7fd882bb-d68f-47b5-8f41-8da40a0b59be
begin
    τ_short_fs = τ_total_frk / 4
    τ_long_fs = τ_total_frk * 2
    feed_fs = Stream(1.0, (I=I0_frk, M=M0_frk))
    branch1_fs, branch2_fs = split_stream(feed_fs, frac_branch1)
    out1_fs = pfr_free_radical(τ_short_fs, k_p_frk, k_d_frk, k_t_frk, f_frk, branch1_fs)
    out2_fs = pfr_free_radical(τ_long_fs, k_p_frk, k_d_frk, k_t_frk, f_frk, branch2_fs)
    combined_fs = mix(out1_fs, out2_fs)
    conv1_fs = 1 - out1_fs.x.M / M0_frk
    conv2_fs = 1 - out2_fs.x.M / M0_frk
    conv_combined_fs = 1 - combined_fs.x.M / M0_frk
end

# ╔═╡ 29e2e8ac-9bbf-4ec5-abb2-6e3ca58d611e
md"""
$(round(100frac_branch1, digits=0))% of the feed to the short-τ branch (τ = $(round(τ_short_fs, sigdigits=3)) s, conversion = $(round(conv1_fs, digits=4))), remainder to the long-τ branch (τ = $(round(τ_long_fs, sigdigits=3)) s, conversion = $(round(conv2_fs, digits=4))).

Remixed product conversion = $(round(conv_combined_fs, digits=4)) — exactly the flow-weighted average of the two branch conversions (a mixer's own definition, `mix`), so it always lies between them.
"""

# ╔═╡ 58deff6b-a56e-4937-9451-149ec44b8324
begin
    fracs_fs = range(0.0, 1.0; length=100)
    convs_combined_fs = Float64[]
    for fr in fracs_fs
        b1, b2 = split_stream(feed_fs, fr)
        o1 = pfr_free_radical(τ_short_fs, k_p_frk, k_d_frk, k_t_frk, f_frk, b1)
        o2 = pfr_free_radical(τ_long_fs, k_p_frk, k_d_frk, k_t_frk, f_frk, b2)
        c = mix(o1, o2)
        push!(convs_combined_fs, 1 - c.x.M / M0_frk)
    end
    plot(fracs_fs, convs_combined_fs;
        xlabel="fraction of feed to short-τ branch", ylabel="remixed conversion",
        label=nothing, lw=2, legend=:right,
        title="Branch-and-remix: conversion vs. split fraction")
    hline!([conv1_fs]; label="short-τ branch alone", ls=:dash, color=:black)
    hline!([conv2_fs]; label="long-τ branch alone", ls=:dot, color=:red)
end

# ╔═╡ fbf992c9-c533-4314-821e-840006f8a295
md"""
## 13. Sanchez-Lacombe mixture activities (consistent, constant-hole-volume)

Section 10 covered mixture *density*; this section covers mixture
*chemical potentials*, via `sl_mixture_activities` (`sanchez_lacombe_activity.jl`),
following von Konigslow, Park & Thompson (2017): the Sanchez-Lacombe
mixture free energy gives thermodynamically consistent chemical
potentials only when the hole volume is held constant with respect to
composition under the derivative — implemented here with its own,
separately verified mixing rules (**not** `sl_mixing_rules` from section
10 — see the module docstring for why those can't be reused). Each
component's activity is reported relative to its own pure fluid at the
same T, P, reusing the same polymer/solvent pair, temperature, and
pressure as sections 1 and 4.
"""

# ╔═╡ 470b50f4-15aa-4118-a112-6de7ea134eea
@bind w1_act Slider(0.02:0.02:0.98; default=0.5, show_value=true)

# ╔═╡ 1587dedf-6585-4f70-96a6-f619a4300feb
begin
    ln_a_act = sl_mixture_ln_activities(T, P_MPa, solv, w1_act, poly, 1 - w1_act)
    a_act = sl_mixture_activities(T, P_MPa, solv, w1_act, poly, 1 - w1_act)
end

# ╔═╡ 3fb2c1ea-02bb-4b80-bc58-3ccb3a57c1c3
md"""
$(round(100w1_act, digits=0))% $(solv.name) / $(round(100 * (1 - w1_act), digits=0))% $(poly.name) by mass, at $(T_C) °C, P = $(P_MPa) MPa:

Solvent ($(solv.name)) activity: ln(a₁) = $(round(ln_a_act.ln_a1, digits=4)), a₁ = $(round(a_act.a1, digits=4))

Polymer ($(poly.name)) activity: ln(a₂) = $(round(ln_a_act.ln_a2, sigdigits=4)) (per mole of chain — astronomically small for a long polymer chain, same "per mole of chain" convention as `ln_activity_polymer` in the Flory-Huggins section)
"""

# ╔═╡ 6f858870-b45a-467f-a6ed-e21b3e88c5dc
begin
    w1s_act = range(0.03, 0.97; length=100)
    ln_a1s_act = [sl_mixture_ln_activities(T, P_MPa, solv, w, poly, 1 - w).ln_a1 for w in w1s_act]
    plot(w1s_act, ln_a1s_act;
        xlabel="w₁ (solvent weight fraction)", ylabel="ln(a₁) (solvent)",
        label=nothing, lw=2,
        title="Sanchez-Lacombe solvent activity vs. composition")
    hline!([0.0]; label="pure solvent reference (ln a₁ = 0)", ls=:dash, color=:black)
end

# ╔═╡ 7001799b-5fa9-4945-bce4-f7c4765aa12f
md"""
## 14. Industrial polyolefin processes: gas-phase, slurry-loop, solution

Polyethylene and polypropylene (`species("polyethylene")`,
`species("polyethylene_ldpe")`, `species("polypropylene")` — section 1)
are made industrially by several reactor configurations that need no new
kinetics beyond `kinetics_coordination.jl`, just the right reactor model
and, for the gas-phase case, converting a monomer *partial pressure* to
the concentration the kinetics expect (`ideal_gas_concentration`). All
three below share the same illustrative propagation/deactivation rate
constants and active-site concentration — only the residence time and
(for the gas-phase case) the monomer supply differ, matching how these
processes actually differ operationally. See the README's process-mapping
table for typical temperature/pressure context (not modeled explicitly
here — this package's coordination kinetics take `k_p`, `k_t` as given,
without a built-in Arrhenius temperature dependence).
"""

# ╔═╡ 59c4baf0-506d-4403-b34d-55815951ad2e
@bind P_gas_MPa Slider(0.5:0.1:3.0; default=1.5, show_value=true)

# ╔═╡ c881a916-4d75-4208-826b-66530076aa50
@bind τ_gas_min Slider(30:10:240; default=120, show_value=true)

# ╔═╡ b4688087-532f-4526-a198-86de57e29143
@bind τ_loop_min Slider(20:5:120; default=75, show_value=true)

# ╔═╡ c3e743fb-7fff-4e34-856f-78c4207b1590
@bind τ_soln_min Slider(0.5:0.5:15; default=3, show_value=true)

# ╔═╡ 6d67a6e6-774e-44e3-a174-922e486d3db2
begin
    k_p_ind = 50.0      # L/(mol s), illustrative coordination propagation rate
    k_t_ind = 5.0e-4    # 1/s, illustrative active-site deactivation rate
    C0_star_ind = 1.0e-4  # mol/L, illustrative active-site concentration
    T_gas_K = 363.15    # 90 °C, illustrative gas-phase temperature (fixed; see note above)

    M_gas_ind = ideal_gas_concentration(P_gas_MPa, T_gas_K)
    res_gas_ind = cstr_coordination(60τ_gas_min, k_p_ind, k_t_ind, C0_star_ind, M_gas_ind)

    M_loop_ind = 3.0    # mol/L, illustrative liquid-phase (diluent) monomer concentration
    res_loop_ind = cstr_coordination(60τ_loop_min, k_p_ind, k_t_ind, C0_star_ind, M_loop_ind)

    M_soln_ind = 2.0    # mol/L, illustrative liquid-phase (solvent) monomer concentration
    res_soln_ind = cstr_coordination(60τ_soln_min, k_p_ind, k_t_ind, C0_star_ind, M_soln_ind)
end

# ╔═╡ 23ef64fb-86c5-4a2d-a7bb-bc45d8680dc6
md"""
Gas-phase: [M] = $(round(M_gas_ind, digits=3)) mol/L from P = $(P_gas_MPa) MPa at $(round(T_gas_K - 273.15, digits=0)) °C, τ = $(τ_gas_min) min → conversion = $(round(res_gas_ind.conversion, digits=4))

Slurry loop: [M] = $(M_loop_ind) mol/L, τ = $(τ_loop_min) min → conversion = $(round(res_loop_ind.conversion, digits=4))

Solution: [M] = $(M_soln_ind) mol/L, τ = $(τ_soln_min) min → conversion = $(round(res_soln_ind.conversion, digits=4))
"""

# ╔═╡ 862b0b7c-9a2f-4536-bc92-5fc30224b44b
bar(["gas-phase", "slurry loop", "solution"],
    [res_gas_ind.conversion, res_loop_ind.conversion, res_soln_ind.conversion];
    ylabel="conversion", legend=nothing,
    title="Conversion by process type (illustrative, same catalyst)")

# ╔═╡ 2b581d13-7845-462f-a41f-04e0f0d73a12
md"""
## 15. LDPE high-pressure autoclave & tubular reactors

The two high-pressure (150-300 MPa) free-radical LDPE processes, both
needing *staged* initiator injection to sustain the radical population
along the reactor (`reactor_staged_injection.jl`): a multi-zone
**autoclave** (each zone an ideal CSTR,
[`cstr_train_free_radical_staged`](@ref)) and a **tubular** reactor with
injection points along its length (each segment a PFR,
[`pfr_train_free_radical_staged`](@ref)). `n` zones/segments share the
same total residence time; `I_first` is charged to the first zone (as
much as could safely be charged in a single dose — dumping the *whole*
charge at once risks a runaway exotherm this package doesn't model), and
an *additional* `I_extra_total` is split evenly across the remaining
zones. Reuses illustrative rate constants of the same order of magnitude
as section 5 (real high-pressure LDPE rate constants differ — see the
note in section 14).

Note: staging the *same* total initiator differently doesn't always
increase conversion (which portion reacts fastest depends on the
nonlinear radical-concentration/decay coupling) — but adding *more* total
initiator, staged in on top of the first dose, never decreases it. That's
the comparison below, and it's also the more industrially accurate one:
staging is what lets a plant add more initiator overall than a single
safe dose would allow, not a free way to redistribute a fixed charge.
"""

# ╔═╡ aecc0a1b-4a01-4786-ad66-3f6a0b09738f
@bind n_zones_ldpe Slider(2:1:6; default=3, show_value=true)

# ╔═╡ 3ee36363-12c5-4212-9249-652823641eaf
@bind I_extra_total_ldpe Slider(0.0:0.002:0.02; default=0.01, show_value=true)

# ╔═╡ 6a0eacdf-1e9a-4028-9f7d-a46700386cd5
begin
    k_p_ldpe, k_d_ldpe, k_t_ldpe, f_ldpe = 1.0e3, 1.0e-4, 1.0e7, 0.5
    I_first_ldpe = 0.01    # mol/L, illustrative first-zone (safe single-dose) charge
    M0_ldpe = 5.0          # mol/L, illustrative bulk-ish monomer concentration
    τ_total_ldpe = 180.0   # s, illustrative total residence time (a few minutes)

    τs_ldpe = fill(τ_total_ldpe / n_zones_ldpe, n_zones_ldpe)
    ΔI_staged_ldpe = vcat(
        [I_first_ldpe],
        fill(I_extra_total_ldpe / (n_zones_ldpe - 1), n_zones_ldpe - 1),
    )
    ΔI_baseline_ldpe = vcat([I_first_ldpe], zeros(n_zones_ldpe - 1))

    res_autoclave = cstr_train_free_radical_staged(τs_ldpe, k_p_ldpe, k_d_ldpe, k_t_ldpe, f_ldpe, ΔI_staged_ldpe, M0_ldpe)
    res_autoclave_baseline = cstr_train_free_radical_staged(τs_ldpe, k_p_ldpe, k_d_ldpe, k_t_ldpe, f_ldpe, ΔI_baseline_ldpe, M0_ldpe)
    res_tubular = pfr_train_free_radical_staged(τs_ldpe, k_p_ldpe, k_d_ldpe, k_t_ldpe, f_ldpe, ΔI_staged_ldpe, M0_ldpe)
    res_tubular_baseline = pfr_train_free_radical_staged(τs_ldpe, k_p_ldpe, k_d_ldpe, k_t_ldpe, f_ldpe, ΔI_baseline_ldpe, M0_ldpe)
end

# ╔═╡ 7ab0c4f0-47a3-4afd-85ad-cd650ae217d0
md"""
With $(n_zones_ldpe) zones/segments, a first-zone charge of $(I_first_ldpe) mol/L, and $(I_extra_total_ldpe) mol/L more split evenly across the rest:

Autoclave: conversion = $(round(res_autoclave.conversion, digits=4)) with staged re-injection vs. $(round(res_autoclave_baseline.conversion, digits=4)) with only the first-zone charge

Tubular: conversion = $(round(res_tubular.conversion, digits=4)) with staged re-injection vs. $(round(res_tubular_baseline.conversion, digits=4)) with only the first-zone charge
"""

# ╔═╡ 37d2758c-8669-447e-8d28-15df9cff47be
begin
    ns_ldpe = 2:6
    convs_staged_auto = Float64[]
    convs_baseline_auto = Float64[]
    for n in ns_ldpe
        τs_n = fill(τ_total_ldpe / n, n)
        ΔI_s = vcat([I_first_ldpe], fill(I_extra_total_ldpe / (n - 1), n - 1))
        ΔI_b = vcat([I_first_ldpe], zeros(n - 1))
        push!(convs_staged_auto, cstr_train_free_radical_staged(τs_n, k_p_ldpe, k_d_ldpe, k_t_ldpe, f_ldpe, ΔI_s, M0_ldpe).conversion)
        push!(convs_baseline_auto, cstr_train_free_radical_staged(τs_n, k_p_ldpe, k_d_ldpe, k_t_ldpe, f_ldpe, ΔI_b, M0_ldpe).conversion)
    end
    plot(ns_ldpe, convs_staged_auto;
        xlabel="number of zones", ylabel="conversion", legend=:bottomright,
        label="first-zone charge + staged re-injection", lw=2, marker=:circle,
        title="Autoclave: staged re-injection vs. first-zone charge alone")
    plot!(ns_ldpe, convs_baseline_auto; label="first-zone charge only", lw=2, ls=:dash, marker=:circle)
end

# ╔═╡ 9523ba48-5bfd-4fca-b0ba-55fa001eaf43
md"""
## 16. Case study: Spheripol-style PP (loop reactors → gas-phase impact copolymer)

A worked example combining several pieces above into one realistic
industrial flowsheet: **two liquid-propylene loop reactors in series**
(the homopolymer/random-copolymer matrix — coordination catalysis, both
loops sharing the same catalyst active sites via
[`cstr_train_coordination`](@ref), since it's the same liquid-propylene
process fluid throughout), followed by **one or two gas-phase reactors**
building an ethylene-propylene rubber (EPR) phase for an *impact*
copolymer grade — the catalyst's active sites (`C_star`) carry over
physically from the loops into the gas phase, but the gas phase is a
genuinely different reactor with its own ethylene/propylene gas feed
(via [`ideal_gas_concentration`](@ref)), **not** a continuation of the
loops' liquid monomer concentration — an earlier draft of this example
got that wrong by chaining `cstr_train_coordination` straight across the
phase change, caught by checking the numbers rather than assumed.
[`instantaneous_copolymer_composition`](@ref) (Mayo-Lewis) gives the
resulting rubber phase's ethylene content from the gas feed composition
and illustrative reactivity ratios — ethylene reacts much faster than
propylene with Ziegler-Natta catalysts, so the rubber ends up markedly
richer in ethylene than the gas feed, the qualitative behavior real
impact-PP processes rely on.
"""

# ╔═╡ a18cd6c9-641e-4fe3-84bb-ca35c6e34a3b
@bind τ_loop_min Slider(20:5:90; default=60, show_value=true)

# ╔═╡ f28ace66-0755-488d-9cd8-62a94e1f2fcc
@bind n_gas_stages Slider(1:1:2; default=1, show_value=true)

# ╔═╡ 1867cf2d-163f-4593-ab41-de9ecddf66e1
@bind P_ethylene_gas Slider(0.1:0.05:0.6; default=0.3, show_value=true)

# ╔═╡ 62d2f338-3689-4257-a9aa-c4d601122988
@bind P_propylene_gas Slider(0.5:0.1:2.0; default=1.2, show_value=true)

# ╔═╡ 7e3ef314-a2fb-442a-8226-5891f3ae49c9
begin
    k_p_loop_ex, k_t_loop_ex = 5.0, 1.0e-4  # illustrative
    C0_star_loop_ex, M_loop_ex = 2.0e-4, 8.0  # mol/L; bulk liquid propylene is concentrated
    τ_loop_s_ex = 60τ_loop_min

    loop_train_ex = cstr_train_coordination(2, τ_loop_s_ex, k_p_loop_ex, k_t_loop_ex, C0_star_loop_ex, M_loop_ex)

    T_gas_ex = 343.15  # K, illustrative gas-phase temperature
    M_ethylene_gas_ex = ideal_gas_concentration(P_ethylene_gas, T_gas_ex)
    M_propylene_gas_ex = ideal_gas_concentration(P_propylene_gas, T_gas_ex)
    M_gas_total_ex = M_ethylene_gas_ex + M_propylene_gas_ex
    τ_gas_s_ex = 3600.0  # s, illustrative per gas-phase stage

    gas_train_ex = cstr_train_coordination(n_gas_stages, τ_gas_s_ex, k_p_loop_ex, k_t_loop_ex, loop_train_ex.C_star, M_gas_total_ex)

    r_ethylene_ex, r_propylene_ex = 5.0, 0.2  # illustrative Ziegler-Natta EP reactivity ratios
    f_ethylene_feed_ex = M_ethylene_gas_ex / M_gas_total_ex
    F_ethylene_rubber_ex = instantaneous_copolymer_composition(r_ethylene_ex, r_propylene_ex, f_ethylene_feed_ex)

    mass_loops_ex = loop_train_ex.conversion * M_loop_ex
    mass_gas_ex = gas_train_ex.conversion * M_gas_total_ex
    rubber_content_ex = mass_gas_ex / (mass_loops_ex + mass_gas_ex)
end

# ╔═╡ 9d84e2a9-4451-4d6c-a5e2-364718e691c0
md"""
**Loops** (2 × $(τ_loop_min) min, propylene): conversion = $(round(loop_train_ex.conversion, digits=4)), active sites remaining = $(round(100loop_train_ex.C_star / C0_star_loop_ex, digits=1))% of the fresh charge

**Gas phase** ($(n_gas_stages) stage(s), C₂⁼ partial pressure $(P_ethylene_gas) MPa, C₃⁼ $(P_propylene_gas) MPa): conversion = $(round(gas_train_ex.conversion, digits=4)), ethylene feed fraction = $(round(100f_ethylene_feed_ex, digits=1))%, ethylene content of the rubber phase = $(round(100F_ethylene_rubber_ex, digits=1))%

**Overall product**: rubber (EPR) content = $(round(100rubber_content_ex, digits=1))% of total polymer produced
"""

# ╔═╡ 98eb5e16-a5bc-44ba-a9eb-49dd54575fe4
begin
    P_ethylenes_ex = range(0.05, 0.6; length=60)
    F_ethylenes_ex = [
        instantaneous_copolymer_composition(
            r_ethylene_ex, r_propylene_ex,
            ideal_gas_concentration(Pe, T_gas_ex) / (ideal_gas_concentration(Pe, T_gas_ex) + M_propylene_gas_ex),
        ) for Pe in P_ethylenes_ex
    ]
    plot(P_ethylenes_ex, F_ethylenes_ex;
        xlabel="ethylene partial pressure (MPa), propylene fixed at $(P_propylene_gas) MPa",
        ylabel="ethylene mole fraction in rubber phase",
        label=nothing, lw=2,
        title="Gas-phase rubber (EPR) composition vs. ethylene feed")
    plot!([0, 0.6], [0, 0.6]; label="rubber = feed (no enrichment)", ls=:dash, color=:black)
end

# ╔═╡ 6f74969a-7a26-435c-ac27-cdefc0037e20
md"""
---
Everything originally on this package's roadmap (recycle loops,
coordination reactors, a general flowsheet solver, and Sanchez-Lacombe
mixture chemical potentials/activities) is now implemented above. See the
repository README for details, including the trade-offs each of these
made along the way.
"""

# ╔═╡ Cell order:
# ╠═1d69aa37-9d2a-4af7-a1d3-9c791a983c12
# ╟─ba17fea8-391c-4a6c-9657-ca1a468b3047
# ╟─bf44a61f-3c6f-4212-8320-922d3dd222dc
# ╠═8546f5e2-21e0-4173-bbad-67fe6a4b4c8d
# ╠═d8a30bc1-4969-4055-9d55-1efc5369b5d6
# ╠═e4e85e1c-02b7-4272-ab76-ba988ccda7d3
# ╠═fd4a18c4-1ad2-466e-b889-98b2c8a5ff68
# ╠═7651a9e8-d0c8-435f-8e1d-21dc48be24dd
# ╟─d9924eaa-de0e-43ee-b9af-eeba20ea4ff2
# ╟─68459002-17d6-4071-89f3-997d3f453ce8
# ╠═664d29e8-9bac-414d-91fb-637e798beb72
# ╠═068f474e-d08c-4204-b6e9-dd61fd2147c0
# ╟─1c1f35f0-9954-40f5-86aa-05fc5da71c63
# ╠═715647ac-d886-47bd-ac78-710f2b0ffc41
# ╟─338160d0-296d-46f3-b760-8f4b2190af0f
# ╠═ab4ee0c7-8a80-4431-81b4-30ef03130da6
# ╟─6103bceb-8a74-49ff-8e03-54a0220ad158
# ╠═66b69bd7-af4e-4465-a0dc-6f167b696b18
# ╠═0600f577-7093-4143-8b69-321c0b2fbd31
# ╠═045bea12-4120-4a56-afc7-b386a98d8936
# ╠═10f6793d-7f12-4f7a-8b39-67a8d4914c7f
# ╟─4b222d62-b1d1-4c02-9f3a-85461ce59b1c
# ╠═fe3637be-868c-4ee2-bc79-899656be2deb
# ╠═151fb00f-1d16-47fc-b605-8f5e095d1665
# ╠═0ad80cab-531b-405d-8837-b9074d32fcc4
# ╠═2689f077-50ab-4917-9797-2b3c8bce3551
# ╟─7534ac44-9487-48af-883d-6e8805b60d70
# ╠═5d467096-e27f-4b7e-9d9b-4c33d26e590e
# ╟─846ed143-a1ae-4fff-b383-e5404b9112bb
# ╠═e5d1ac62-d0e6-4efa-aa5a-302ea607b969
# ╠═12428167-5594-4efa-bef9-175a32d67f1d
# ╠═69263b58-861a-41a0-81bd-6718c32ad6fb
# ╠═52849a69-c1ac-418e-8f52-35958b930477
# ╟─f2306b92-42be-4d25-a276-55e45bb8401d
# ╠═38670f3b-dd5a-4a13-91dd-49db4586e034
# ╠═0028e23f-fa21-4be0-a9f9-51a727765eff
# ╟─065d37c5-ab38-430f-a26b-3682ea251778
# ╠═f6c82001-9a3d-49ab-871b-3a8c77f26a95
# ╠═3ae19407-0464-49c4-a3f1-0b305be274d2
# ╟─d9972239-b3e5-4eba-bd55-d4f90c770282
# ╠═7ad99fa2-83f8-450f-8f70-e4dcc2c6edfb
# ╠═48b3199c-d5c5-44a6-8ef6-fd93e59d4e78
# ╠═6f8f3320-8e35-4c43-b7c4-6f5d2d4ff11f
# ╠═1e0ff82a-d49e-4c10-bf22-a0fc4edcd9fb
# ╟─56ba9de1-48e7-4252-9b7d-dc4a1e7900dd
# ╠═da08bb91-9401-4b31-9952-ba1fccf997ff
# ╠═4bd841d5-5e7d-48c8-994a-2f48bc049dc5
# ╟─8139da23-ae9c-4052-8bd3-32ef4fa0391f
# ╠═df12ac0e-1ac5-4ed5-80bc-78be20d8940c
# ╠═fefabb16-4409-48a6-b6d6-10ae48e15e97
# ╟─142995ab-27e7-45fe-9cd4-2e7a6bd60151
# ╠═91358a8b-c845-4c13-9779-d9d01ac95cd6
# ╟─0111a035-58a0-45ed-8913-0fea5774366e
# ╠═bba86433-0285-4076-8b58-58e6c870362d
# ╠═3977b7bf-5ee5-49b0-a831-9a9332e7e73d
# ╟─9497e799-3098-446c-892e-b485ac224e9a
# ╠═6712044c-4a44-4e58-b3a9-d00604247b2f
# ╟─73fbbb0f-2400-43e9-8e9e-4a9161a21414
# ╠═8822e6f6-0f04-420c-9387-70607c8e6cf4
# ╠═569b1802-d2fc-4092-a694-11b6703fc640
# ╟─dd2436b1-9501-4982-9d6b-3f3649c45100
# ╠═66ee83fc-148a-4a85-91a0-fbba76f70197
# ╟─178d60eb-0691-4815-8822-d1d98ce1a893
# ╠═a2a2b30c-5670-48f4-9ba9-54fcd576e200
# ╠═7fd882bb-d68f-47b5-8f41-8da40a0b59be
# ╟─29e2e8ac-9bbf-4ec5-abb2-6e3ca58d611e
# ╠═58deff6b-a56e-4937-9451-149ec44b8324
# ╟─fbf992c9-c533-4314-821e-840006f8a295
# ╠═470b50f4-15aa-4118-a112-6de7ea134eea
# ╠═1587dedf-6585-4f70-96a6-f619a4300feb
# ╟─3fb2c1ea-02bb-4b80-bc58-3ccb3a57c1c3
# ╠═6f858870-b45a-467f-a6ed-e21b3e88c5dc
# ╟─7001799b-5fa9-4945-bce4-f7c4765aa12f
# ╠═59c4baf0-506d-4403-b34d-55815951ad2e
# ╠═c881a916-4d75-4208-826b-66530076aa50
# ╠═b4688087-532f-4526-a198-86de57e29143
# ╠═c3e743fb-7fff-4e34-856f-78c4207b1590
# ╠═6d67a6e6-774e-44e3-a174-922e486d3db2
# ╟─23ef64fb-86c5-4a2d-a7bb-bc45d8680dc6
# ╠═862b0b7c-9a2f-4536-bc92-5fc30224b44b
# ╟─2b581d13-7845-462f-a41f-04e0f0d73a12
# ╠═aecc0a1b-4a01-4786-ad66-3f6a0b09738f
# ╠═3ee36363-12c5-4212-9249-652823641eaf
# ╠═6a0eacdf-1e9a-4028-9f7d-a46700386cd5
# ╟─7ab0c4f0-47a3-4afd-85ad-cd650ae217d0
# ╠═37d2758c-8669-447e-8d28-15df9cff47be
# ╟─9523ba48-5bfd-4fca-b0ba-55fa001eaf43
# ╠═a18cd6c9-641e-4fe3-84bb-ca35c6e34a3b
# ╠═f28ace66-0755-488d-9cd8-62a94e1f2fcc
# ╠═1867cf2d-163f-4593-ab41-de9ecddf66e1
# ╠═62d2f338-3689-4257-a9aa-c4d601122988
# ╠═7e3ef314-a2fb-442a-8226-5891f3ae49c9
# ╟─9d84e2a9-4451-4d6c-a5e2-364718e691c0
# ╠═98eb5e16-a5bc-44ba-a9eb-49dd54575fe4
# ╟─6f74969a-7a26-435c-ac27-cdefc0037e20
