from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

from symplyphysics.definitions import volume_number_density

# Description
## The atomic number density (N; atoms/cm^3) is the number of atoms of a given type per unit volume (V; cm^3)
## of the material.

## Law: N = ⍴ * Na / M
## Where:
## ⍴ (material density) is density of the material
## Na is Avogadro's number
## M (atomic or molecular weight) - total weight of an atom
## N is the atomic number density

atomic_number_density, material_density, atomic_weight = symbols('atomic_number_density material_density atomic_weight')
law = Eq(atomic_number_density, material_density * units.avogadro / atomic_weight)

# Derive the same law from volume number density law
material_mass = symbols('material_mass')

#TODO: import from definitions
density_law = Eq(material_density, material_mass / volume_number_density.volume)
avogadro_law = Eq(volume_number_density.objects, material_mass * units.avogadro / atomic_weight)

derived_law = [volume_number_density.definition, density_law, avogadro_law]

## Check the equivalence of 'law' and 'derived_law'
derived_number_density = solve(derived_law,
    (material_mass, volume_number_density.objects, volume_number_density.number_density),
    dict=True)[0][volume_number_density.number_density]
assert solve(law, atomic_number_density, dict=True)[0][atomic_number_density] == derived_number_density

def print():
    return pretty(law, use_unicode=False)

@validate_input(material_density_=(units.mass / units.volume), atomic_weight_=(units.mass / units.amount_of_substance))
@validate_output(1 / units.volume)
def calculate_atomic_number_density(material_density_: Quantity, atomic_weight_: Quantity) -> Quantity:
    solved = solve(law, atomic_number_density, dict=True)[0][atomic_number_density]
    result_expr = solved.subs({
        material_density: material_density_,
        atomic_weight: atomic_weight_})
    return expr_to_quantity(result_expr, 'atomic_number_density')