r"""
Number density via volumetric density and molar mass
====================================================

Number density, or number of molecules per unit volume, can be expressed using
the volumetric density, or mass per unit volume, and molar mass, or mass per unit
amount of substance.

**Notation:**

#. :quantity_notation:`avogadro_constant`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Number_density#Mass_density>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, quantities, symbols
from symplyphysics.definitions import number_density_is_number_of_objects_per_unit_volume
from symplyphysics.definitions import density_from_mass_volume
from symplyphysics.laws.chemistry import avogadro_constant_is_particle_count_over_amount_of_substance
from symplyphysics.laws.quantities import quantity_is_molar_quantity_times_amount_of_substance as molar_qty_law

number_density = symbols.number_density
"""
:symbols:`number_density`, or number of molecules per unit volume.
"""

volumetric_density = symbols.density
r"""
Volumetric :symbols:`density`, or mass per unit volume.
"""

molar_mass = symbols.molar_mass
"""
:symbols:`molar_mass`, or mass per unit amount of substance.
"""

law = Eq(number_density, volumetric_density * quantities.avogadro_constant / molar_mass)
"""
:laws:symbol::

:laws:latex::
"""

# Derive the same law from volume number density law

_density_law = density_from_mass_volume.definition.subs({
    density_from_mass_volume.volume: number_density_is_number_of_objects_per_unit_volume.volume,
    density_from_mass_volume.density: volumetric_density
})

_avogadro_law = avogadro_constant_is_particle_count_over_amount_of_substance.law.subs({
    avogadro_constant_is_particle_count_over_amount_of_substance.particle_count:
    number_density_is_number_of_objects_per_unit_volume.number_of_objects
})

_atomic_weight_law = molar_qty_law.law.subs({
    molar_qty_law.molar_quantity:
        molar_mass,
    molar_qty_law.extensive_quantity:
    density_from_mass_volume.mass,
    molar_qty_law.amount_of_substance:
    avogadro_constant_is_particle_count_over_amount_of_substance.amount_of_substance
})

_derived_law = [
    number_density_is_number_of_objects_per_unit_volume.definition, _density_law, _avogadro_law,
    _atomic_weight_law
]

## Check the equivalence of 'law' and '_derived_law'
_derived_number_density = solve(_derived_law, (density_from_mass_volume.mass,
    number_density_is_number_of_objects_per_unit_volume.number_of_objects,
    number_density_is_number_of_objects_per_unit_volume.number_density,
    avogadro_constant_is_particle_count_over_amount_of_substance.amount_of_substance),
    dict=True)[0][number_density_is_number_of_objects_per_unit_volume.number_density]
assert solve(law, number_density, dict=True)[0][number_density] == _derived_number_density


@validate_input(material_density_=volumetric_density, atomic_weight_=molar_mass)
@validate_output(number_density)
def calculate_atomic_number_density(material_density_: Quantity,
    atomic_weight_: Quantity) -> Quantity:
    solved = solve(law, number_density, dict=True)[0][number_density]
    result_expr = solved.subs({volumetric_density: material_density_, molar_mass: atomic_weight_})
    return Quantity(result_expr)
