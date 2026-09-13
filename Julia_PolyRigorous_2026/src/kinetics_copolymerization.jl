"""
Copolymer composition: the Mayo-Lewis equation for the instantaneous
composition of a copolymer chain from the comonomer feed composition and
the two reactivity ratios — standard textbook material (Odian,
*Principles of Polymerization*, Ch. 6), implemented directly.

Convention: monomers are labeled 1 and 2; `f1 = [M1]/([M1]+[M2])` is their
feed mole fraction, and `r1 = k11/k12`, `r2 = k22/k21` the reactivity
ratios (`kij` the rate constant for a chain ending in monomer `i` adding
monomer `j`) — `r1 > 1` means a chain ending in monomer 1 prefers adding
another monomer 1 over monomer 2, and vice versa.

Scope note: this covers *instantaneous composition* only — not
composition drift with conversion (which needs integrating the
Mayo-Lewis equation as the feed depletes and is composition-dependent
on which monomer is being consumed faster) and not copolymerization
*rate* (which needs the individual `kij`, not just their ratios). For a
reactor run continuously at (approximately) steady feed composition —
the relevant case for the CSTR-based industrial examples in this
package (see [`cstr_coordination`](@ref) and the Spheripol notebook
example) — the instantaneous composition *is* the steady-state product
composition, since the feed composition seen by the catalyst doesn't
drift the way it would in a batch reactor.
"""

"""
    instantaneous_copolymer_composition(r1, r2, f1)

`F1`, the instantaneous mole fraction of monomer 1 incorporated into the
copolymer, from the Mayo-Lewis equation:

``F_1 = \\dfrac{r_1 f_1^2 + f_1 f_2}{r_1 f_1^2 + 2 f_1 f_2 + r_2 f_2^2}, \\qquad f_2 = 1 - f_1``

Reduces to `F1 = f1` (the copolymer composition matches the feed exactly,
"ideal"/Bernoullian copolymerization) when `r1 = r2 = 1`, and to pure
monomer-1 (`F1 = 1`) / monomer-2 (`F1 = 0`) incorporation at the
corresponding feed limits — both checked in the test suite, along with
the identity `F1(r1, r2, f1) + F1(r2, r1, 1 - f1) = 1` (monomer 2's own
instantaneous fraction is just `1 - F1`, computed the same way with the
monomers' roles swapped).
"""
function instantaneous_copolymer_composition(r1, r2, f1)
    0 <= f1 <= 1 || throw(ArgumentError("f1 must be in [0, 1]"))
    f2 = 1 - f1
    return (r1 * f1^2 + f1 * f2) / (r1 * f1^2 + 2 * f1 * f2 + r2 * f2^2)
end

"""
    azeotrope_composition(r1, r2)

The feed mole fraction `f1` at which the instantaneous copolymer
composition equals the feed composition (`F1 = f1`), so it doesn't drift
with conversion even in a batch/depleting-feed reactor:

``f_{1,azeo} = \\dfrac{1 - r_2}{2 - r_1 - r_2}``

Only meaningful (in `(0, 1)`) when `r1` and `r2` are both less than 1 or
both greater than 1; throws for `r1 + r2 = 2` (no azeotrope, or every
composition is one, as in the ideal `r1 = r2 = 1` case) rather than
silently returning a meaningless value.
"""
function azeotrope_composition(r1, r2)
    denom = 2 - r1 - r2
    denom != 0 || throw(ArgumentError("no well-defined azeotrope when r1 + r2 = 2"))
    return (1 - r2) / denom
end
