"""
Enthalpy of dielectrics
=======================

The enthalpy of a system with a dielectric medium can be expressed using
its internal energy and the electric field and electric displacement.

**Conditions:**

#. The dielectric is isotropic whether or not the electric field is present.
#. The medium is homogeneous.
#. The volume change of the medium is insignificant.

**Links:**

#. Formula 31.5 on p. 122 of "General Course of Physics" (Obschiy kurs fiziki), vol. 3 by Sivukhin D.V. (1979).
"""

from sympy import Eq
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

enthalpy_density = Symbol("H", units.energy / units.volume)
"""
:symbols:`enthalpy` of the system per unit :symbols:`volume`.
"""

internal_energy_density = Symbol("U", units.energy / units.volume)
"""
:symbols:`internal_energy` of the system per units :symbols:`volume`.
"""

electric_field_strength = symbols.electric_field_strength
"""
:symbols:`electric_field_strength` in the medium.
"""

electric_displacement = symbols.electric_displacement
"""
:symbols:`electric_displacement` in the medium.
"""

law = Eq(
    enthalpy_density,
    internal_energy_density - electric_field_strength * electric_displacement,
)


@validate_input(
    internal_energy_density_=internal_energy_density,
    electric_field_strength_=electric_field_strength,
    electric_displacement_=electric_displacement,
)
@validate_output(enthalpy_density)
def calculate_enthalpy_density(
    internal_energy_density_: Quantity,
    electric_field_strength_: Quantity,
    electric_displacement_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        internal_energy_density: internal_energy_density_,
        electric_field_strength: electric_field_strength_,
        electric_displacement: electric_displacement_,
    })
    return Quantity(result)
