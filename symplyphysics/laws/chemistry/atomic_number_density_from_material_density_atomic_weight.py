from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)
from symplyphysics.definitions import volume_number_density
from symplyphysics.definitions import density_from_mass_volume
from symplyphysics.laws.chemistry import avogadro_number_from_mole_count
from symplyphysics.laws.chemistry import atomic_weight_from_mass_mole_count

# Description
## The atomic number density (N; atoms/cm^3) is the number of atoms of a given type per unit volume (V; cm^3)
## of the material.

## Law: N = ⍴ * Na / M
## Where:
## ⍴ (material density) is density of the material.
##   See [material density](symplyphysics/definitions/density_from_mass_volume.py) implementation.
## Na is Avogadro's number.
##   See [avogadro number](./avogadro_number_from_mole_count.py) implementation.
## M (atomic or molecular weight) - total weight of an atom.
##   See [atomic weight](./atomic_weight_from_mass_mole_count.py) implementation.
## N is the atomic number density.

atomic_number_density = Symbol("atomic_number_density", 1 / units.volume)
material_density = Symbol("material_density", units.mass / units.volume)
atomic_weight = Symbol("atomic_weight", units.mass / units.amount_of_substance)

law = Eq(atomic_number_density, material_density * units.avogadro / atomic_weight)

# Derive the same law from volume number density law

density_law = density_from_mass_volume.definition.subs({
    density_from_mass_volume.volume: volume_number_density.volume,
    density_from_mass_volume.density: material_density
})

avogadro_law = avogadro_number_from_mole_count.law.subs(
    {avogadro_number_from_mole_count.particles_count: volume_number_density.objects})

atomic_weight_law = atomic_weight_from_mass_mole_count.law.subs({
    atomic_weight_from_mass_mole_count.atomic_weight:
        atomic_weight,
    atomic_weight_from_mass_mole_count.symbols.basic.mass:
    density_from_mass_volume.symbols.basic.mass,
    atomic_weight_from_mass_mole_count.mole_count:
    avogadro_number_from_mole_count.mole_count
})

derived_law = [volume_number_density.definition, density_law, avogadro_law, atomic_weight_law]

## Check the equivalence of 'law' and 'derived_law'
derived_number_density = solve(derived_law,
    (density_from_mass_volume.symbols.basic.mass, volume_number_density.objects,
    volume_number_density.number_density, avogadro_number_from_mole_count.mole_count),
    dict=True)[0][volume_number_density.number_density]
assert solve(law, atomic_number_density,
    dict=True)[0][atomic_number_density] == derived_number_density


def print_law() -> str:
    return print_expression(law)


@validate_input(material_density_=material_density, atomic_weight_=atomic_weight)
@validate_output(atomic_number_density)
def calculate_atomic_number_density(material_density_: Quantity,
    atomic_weight_: Quantity) -> Quantity:
    solved = solve(law, atomic_number_density, dict=True)[0][atomic_number_density]
    result_expr = solved.subs({material_density: material_density_, atomic_weight: atomic_weight_})
    return Quantity(result_expr)
