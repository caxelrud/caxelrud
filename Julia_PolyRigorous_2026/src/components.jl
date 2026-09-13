"""
    Species

A pure-component record holding both classical solution-thermodynamics
parameters (molar mass, molar volume, Hildebrand solubility parameter) and
Sanchez-Lacombe lattice-fluid characteristic parameters (T*, P*, ρ*).

Fields
- `name`      :: component name
- `kind`      :: `:solvent` or `:polymer`
- `M`         :: molar mass, g/mol. For polymers this is the number-average
                 molar mass `Mn` of the *chain used in the calculation*, not
                 an intrinsic constant of the polymer species.
- `Vm`        :: molar volume at the reference conditions, cm^3/mol
- `δ`         :: Hildebrand solubility parameter, MPa^0.5
- `Tstar`     :: Sanchez-Lacombe characteristic temperature, K
- `Pstar`     :: Sanchez-Lacombe characteristic pressure, MPa
- `rhostar`   :: Sanchez-Lacombe characteristic (close-packed) density, kg/m^3
"""
struct Species
    name::String
    kind::Symbol
    M::Float64
    Vm::Float64
    δ::Float64
    Tstar::Float64
    Pstar::Float64
    rhostar::Float64
end

function Species(; name, kind, M, Vm, δ, Tstar, Pstar, rhostar)
    kind in (:solvent, :polymer) || throw(ArgumentError("kind must be :solvent or :polymer"))
    return Species(name, kind, M, Vm, δ, Tstar, Pstar, rhostar)
end

"""
    degree_of_polymerization(polymer::Species, solvent::Species)

Flory-Huggins size ratio `N = Vm(polymer)/Vm(solvent)` between one polymer
chain (of molar mass `polymer.M`, i.e. the `Mn` you constructed it with) and
one solvent molecule. This is what enters the Flory-Huggins expressions as
`n`/`N`.
"""
degree_of_polymerization(polymer::Species, solvent::Species) = polymer.Vm / solvent.Vm

"""
    segment_number(sp::Species)

Sanchez-Lacombe number of lattice sites per molecule,
`r = M * Pstar / (R * Tstar * rhostar)`, with `M` in kg/mol, `Pstar` in Pa,
`rhostar` in kg/m^3 and `R` the gas constant in SI units. Dimensionless.
"""
function segment_number(sp::Species)
    R = 8.314462618       # J/(mol K)
    M_kg = sp.M / 1000     # kg/mol
    Pstar_Pa = sp.Pstar * 1e6
    return M_kg * Pstar_Pa / (R * sp.Tstar * sp.rhostar)
end

# ----------------------------------------------------------------------------
# Built-in illustrative species database.
#
# NOTE ON DATA PROVENANCE: the Sanchez-Lacombe (Tstar, Pstar, rhostar) values
# below are typical literature-compilation figures (of the kind tabulated in
# reviews such as Rodgers, J. Appl. Polym. Sci. 1993, and Sanchez & Panayiotou,
# "Models for Thermodynamic and Phase Equilibria Calculations", 1994) rounded
# to 3-4 significant figures. They are meant to make the package runnable
# out of the box and to be *representative*, not to be treated as validated
# regression parameters for a specific resin grade or solvent lot. For real
# engineering work, refit (Tstar, Pstar, rhostar) against PVT data for your
# actual material before trusting the numbers.
#
# "polyethylene_ldpe" and "polypropylene" specifically use the pure-fluid
# Sanchez-Lacombe parameters (Table 1) from von Konigslow, Thompson et al.,
# *Soft Matter* 14, 4603 (2018) — regressed from real PVT data (LDPE: Hasan
# et al.; linear PP: same experimental program) rather than a generic
# compilation figure, since that paper's PDF was already read directly
# (rendered to images, not text-extracted, to avoid garbling the parameter
# table) while researching the mixture chemical-potential work above.
# `Vm` (an actual, not close-packed, reference density) and `δ` are typical
# commodity-polyolefin literature values, not from that paper (which doesn't
# report them) — polyolefins are chemically similar enough that these don't
# vary much between grades.
# ----------------------------------------------------------------------------

const SPECIES_DB = Dict{String,Species}(
    "polystyrene" => Species(
        name="Polystyrene", kind=:polymer, M=100_000.0, Vm=95_000.0,
        δ=18.6, Tstar=735.0, Pstar=357.0, rhostar=1105.0,
    ),
    "pmma" => Species(
        name="Poly(methyl methacrylate)", kind=:polymer, M=100_000.0, Vm=78_600.0,
        δ=19.4, Tstar=696.0, Pstar=503.0, rhostar=1269.0,
    ),
    "polyethylene" => Species(
        name="Polyethylene (HDPE-like)", kind=:polymer, M=100_000.0, Vm=110_000.0,
        δ=16.9, Tstar=649.0, Pstar=425.0, rhostar=904.0,
    ),
    "polyethylene_ldpe" => Species(
        name="Polyethylene (LDPE)", kind=:polymer, M=100_000.0, Vm=108_900.0,
        δ=16.9, Tstar=586.6, Pstar=407.5, rhostar=927.1,
    ),
    "polypropylene" => Species(
        name="Polypropylene (isotactic, linear)", kind=:polymer, M=100_000.0, Vm=110_500.0,
        δ=16.8, Tstar=662.8, Pstar=316.2, rhostar=868.5,
    ),
    "toluene" => Species(
        name="Toluene", kind=:solvent, M=92.14, Vm=106.9,
        δ=18.2, Tstar=635.0, Pstar=419.0, rhostar=969.0,
    ),
    "benzene" => Species(
        name="Benzene", kind=:solvent, M=78.11, Vm=89.4,
        δ=18.7, Tstar=630.0, Pstar=444.0, rhostar=999.0,
    ),
    "cyclohexane" => Species(
        name="Cyclohexane", kind=:solvent, M=84.16, Vm=108.7,
        δ=16.8, Tstar=650.0, Pstar=354.0, rhostar=896.0,
    ),
    "n-hexane" => Species(
        name="n-Hexane", kind=:solvent, M=86.18, Vm=131.6,
        δ=14.9, Tstar=476.0, Pstar=298.0, rhostar=775.0,
    ),
)

"""
    species(key::AbstractString)

Look up a built-in `Species` by key (case-insensitive), e.g.
`species("polystyrene")`, `species("toluene")`. See `SPECIES_DB` for the
full list of keys.
"""
function species(key::AbstractString)
    k = lowercase(key)
    haskey(SPECIES_DB, k) || throw(KeyError("Unknown species \"$key\". Available: $(join(sort(collect(keys(SPECIES_DB))), ", "))"))
    return SPECIES_DB[k]
end
