"""
Step-growth (condensation) polymerization kinetics and the resulting
Flory "most probable" molecular weight distribution.

Convention follows Odian, *Principles of Polymerization*, Ch. 2: `p` is the
extent of reaction (fraction of functional groups reacted), `c0` the
initial concentration of the limiting functional group, and `k` the
relevant rate constant.
"""

"""
    extent_reaction_external_catalyst(k, c0, t)

Extent of reaction `p(t)` for step-growth polymerization with an external
catalyst held at constant concentration, so the kinetics are second order
overall (first order in each of the two reacting functional groups):
`-dc/dt = k c^2`, integrating to `1/(1-p) = 1 + k c0 t`, i.e.

``p(t) = \\dfrac{k c_0 t}{1 + k c_0 t}``
"""
extent_reaction_external_catalyst(k, c0, t) = (k * c0 * t) / (1 + k * c0 * t)

"""
    extent_reaction_self_catalyzed(k, c0, t)

Extent of reaction `p(t)` for step-growth polymerization with no external
catalyst, where the reacting functional group (e.g. a carboxylic acid)
catalyzes its own reaction, making the kinetics third order overall:
`-dc/dt = k c^3`, integrating to `1/(1-p)^2 = 1 + 2 k c0^2 t`, i.e.

``p(t) = 1 - \\dfrac{1}{\\sqrt{1 + 2 k c_0^2 t}}``
"""
extent_reaction_self_catalyzed(k, c0, t) = 1 - 1 / sqrt(1 + 2 * k * c0^2 * t)

"""
    carothers_Xn(p; r=1.0)

Number-average degree of polymerization from the extent of reaction `p`
(Carothers equation). `r` is the stoichiometric ratio of the two
functional groups (or, with a monofunctional impurity present, the
ratio accounting for it); `r = 1.0` (the default) is the stoichiometrically
balanced case `Xn = 1/(1-p)`. For `r < 1`:

``X_n = \\dfrac{1 + r}{1 + r - 2 r p}``
"""
function carothers_Xn(p; r=1.0)
    0 <= p < 1 || throw(ArgumentError("p must be in [0, 1)"))
    0 < r <= 1 || throw(ArgumentError("r must be in (0, 1]"))
    return (1 + r) / (1 + r - 2 * r * p)
end

"""
    flory_mole_fraction(x, p)

Mole fraction of `x`-mers in the Flory "most probable" distribution:
`n_x = (1-p) p^(x-1)`.
"""
flory_mole_fraction(x, p) = (1 - p) * p^(x - 1)

"""
    flory_weight_fraction(x, p)

Weight fraction of `x`-mers in the Flory "most probable" distribution:
`w_x = x (1-p)^2 p^(x-1)`.
"""
flory_weight_fraction(x, p) = x * (1 - p)^2 * p^(x - 1)

"""
    flory_Xw(p)

Weight-average degree of polymerization of the Flory distribution:
`Xw = (1+p)/(1-p)`.
"""
flory_Xw(p) = (1 + p) / (1 - p)

"""
    flory_PDI(p)

Polydispersity index `Xw/Xn = 1 + p` of the Flory distribution; approaches
2 as `p -> 1` (the classic step-growth result).
"""
flory_PDI(p) = 1 + p
