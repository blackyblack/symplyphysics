"""
Gibbs energy of dielectrics
===========================

Gibbs energy of the system with a dielectric medium can be expressed using
the Helmholtz free energy and the electric displacement and field strength.

**Conditions:**

#. The dielectric is isotropic whether or not the electric field is present.
#. The medium is homogeneous.
#. The volume change of the medium is insignificant.

**Links:**

#. Formula 31.4 on p. 122 of "General Course of Physics" (Obschiy kurs fiziki), vol. 3 by Sivukhin D.V. (1979).
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

gibbs_energy_density = Symbol("G", units.energy / units.volume)
"""
:symbols:`gibbs_energy` of the system per unit :symbols:`volume`.
"""

free_energy_density = Symbol("F", units.energy / units.volume)
"""
:symbols:`helmholtz_free_energy` of the system per units :symbols:`volume`.
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
    gibbs_energy_density,
    free_energy_density - electric_field_strength * electric_displacement,
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    free_energy_density_=free_energy_density,
    electric_field_strength_=electric_field_strength,
    electric_displacement_=electric_displacement,
)
@validate_output(gibbs_energy_density)
def calculate_gibbs_energy_density(
    free_energy_density_: Quantity,
    electric_field_strength_: Quantity,
    electric_displacement_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        free_energy_density: free_energy_density_,
        electric_field_strength: electric_field_strength_,
        electric_displacement: electric_displacement_,
    })
    return Quantity(result)
