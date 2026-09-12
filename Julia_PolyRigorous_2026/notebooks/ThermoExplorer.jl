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
@bind chi_s Slider(0.20:0.01:0.50; default=0.34, show_value=true)

# ╔═╡ 068f474e-d08c-4204-b6e9-dd61fd2147c0
begin
    chi = chi_from_solubility(solv.Vm, solv.delta, poly.delta, T; chi_s=chi_s)
    crit = critical_point(N)
end

# ╔═╡ 1c1f35f0-9954-40f5-86aa-05fc5da71c63
md"""
χ (from solubility parameters, entropic correction χ_s = $(chi_s)) = **$(round(chi, digits=3))**

Critical point for this N: φ₂,c = $(round(crit.phi2c, digits=4)), χ_c = $(round(crit.chic, digits=3))

$(chi > crit.chic ? "⚠️ χ > χ_c at this temperature: the model predicts liquid-liquid phase separation over part of the composition range." : "✅ χ < χ_c at this temperature: the model predicts complete miscibility at all compositions.")
"""

# ╔═╡ 715647ac-d886-47bd-ac78-710f2b0ffc41
begin
    phis = range(0.001, 0.999; length=300)
    ln_a1 = [ln_activity_solvent(p, N, chi) for p in phis]
    plot(phis, ln_a1;
        xlabel="polymer volume fraction φ₂", ylabel="ln(a₁)",
        label="ln(solvent activity)", lw=2, legend=:bottomleft,
        title="Solvent activity vs. composition")
end

# ╔═╡ 338160d0-296d-46f3-b760-8f4b2190af0f
md"## 3. Phase diagram (χ vs. φ₂)"

# ╔═╡ ab4ee0c7-8a80-4431-81b4-30ef03130da6
begin
    spin = spinodal_curve(N; npoints=300)
    plt = plot([s.phi2 for s in spin], [s.chi for s in spin];
        label="spinodal", lw=2, xlabel="φ₂", ylabel="χ",
        title="Phase diagram at N = $(round(N, digits=1))", legend=:topright)
    scatter!(plt, [crit.phi2c], [crit.chic]; label="critical point", ms=6)
    if chi > crit.chic
        bcurve = binodal_curve(N; chi_max=max(1.001chi, 1.2crit.chic), npoints=40)
        if !isempty(bcurve)
            plot!(plt, [b.phi2a for b in bcurve], [b.chi for b in bcurve];
                label="binodal", lw=2, ls=:dash, color=3)
            plot!(plt, [b.phi2b for b in bcurve], [b.chi for b in bcurve];
                label=nothing, lw=2, ls=:dash, color=3)
        end
    end
    hline!(plt, [chi]; label="current χ", ls=:dot, color=:black)
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

# ╔═╡ bf5f98f5-b671-49b3-9f3b-e4b01165994e
md"""
---
**Roadmap** (not yet implemented in this version): polymerization kinetics
(free-radical / step-growth / coordination) for molecular weight
distributions, reactor unit operations (CSTR/PFR/batch) built on those
kinetics, Sanchez-Lacombe *mixture* thermodynamics (binary mixing rules and
chemical potentials), and eventually a full flowsheet solver. See the
repository README for details.
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
# ╟─bf5f98f5-b671-49b3-9f3b-e4b01165994e
