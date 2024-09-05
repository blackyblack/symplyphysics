r"""
Number density via volumetric density and molar mass
====================================================

Number density, or number of molecules per unit volume, can be expressed using
the volumetric density, or mass per unit volume, and molar mass, or mass per unit
amount or substance.

**Notation:**

#. :math:`N_\text{A}` (:code:`N_A`) is the Avogadro constant.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)
from symplyphysics.definitions import number_density_is_number_of_objects_per_unit_volume
from symplyphysics.definitions import density_from_mass_volume
from symplyphysics.laws.chemistry import avogadro_constant_is_particle_count_over_amount_of_substance
from symplyphysics.laws.quantities import quantity_is_molar_quantity_times_amount_of_substance as molar_qty_law

number_density = Symbol("number_density", 1 / units.volume)
"""
Number density, or number of molecules per unit volume.

Symbol:
    :code:`n`
"""

volumetric_density = Symbol("volumetric_density", units.mass / units.volume)
r"""
Volumetric density, or mass per unit volume.

Symbol:
    :code:`rho`

Latex:
    :math:`\rho`
"""

molar_mass = Symbol("molar_mass", units.mass / units.amount_of_substance)
"""
Molar mass, or mass per unit amount of substance.

Symbol:
    :code:`M`
"""

law = Eq(number_density, volumetric_density * units.avogadro / molar_mass)
r"""
:code:`n = rho * N_A / M`

Latex:
    .. math::
        n = \frac{\rho N_\text{A}}{M}
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
    molar_qty_law.molar_quantity: molar_mass,
    molar_qty_law.extensive_quantity: density_from_mass_volume.mass,
    molar_qty_law.amount_of_substance: avogadro_constant_is_particle_count_over_amount_of_substance.amount_of_substance
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
assert solve(law, number_density,
    dict=True)[0][number_density] == _derived_number_density


@validate_input(material_density_=volumetric_density, atomic_weight_=molar_mass)
@validate_output(number_density)
def calculate_atomic_number_density(material_density_: Quantity,
    atomic_weight_: Quantity) -> Quantity:
    solved = solve(law, number_density, dict=True)[0][number_density]
    result_expr = solved.subs({volumetric_density: material_density_, molar_mass: atomic_weight_})
    return Quantity(result_expr)
