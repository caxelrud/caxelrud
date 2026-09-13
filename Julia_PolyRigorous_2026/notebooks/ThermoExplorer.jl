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

# ╔═╡ 6f74969a-7a26-435c-ac27-cdefc0037e20
md"""
---
**Roadmap** (not yet implemented in this version): Sanchez-Lacombe
*mixture* thermodynamics (binary mixing rules and chemical potentials),
recycle loops, and eventually a full flowsheet solver. See the repository
README for details.
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
# ╟─6f74969a-7a26-435c-ac27-cdefc0037e20
