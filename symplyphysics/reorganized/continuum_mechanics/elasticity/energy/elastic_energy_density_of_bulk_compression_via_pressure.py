"""
Elastic energy density of bulk compression via pressure
=======================================================

Volumetric density of the elastic energy of a body that is under bulk
compression is proportional to the square of the pressure in the body
and inversely proportional to the bulk modulus of the body's material.

**Conditions:**

#. The body undergoes bulk compression.

**Links:**

#. Formula 77.5 on p. 393 of "General Course of Physics" (Obschiy kurs fiziki), vol. 1 by Sivukhin D.V. (1979).
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

elastic_energy_density = symbols.energy_density
"""
Elastic energy of the deformed body per unit of its volume. See :symbols:`energy_density`
"""

pressure = symbols.pressure
"""
:symbols:`pressure` in the deformed body.
"""

bulk_modulus = symbols.bulk_modulus
"""
:symbols:`bulk_modulus` of the material.
"""

law = Eq(elastic_energy_density, pressure**2 / (2 * bulk_modulus))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    pressure_=pressure,
    bulk_modulus_=bulk_modulus,
)
@validate_output(elastic_energy_density)
def calculate_elastic_energy_density(
    pressure_: Quantity,
    bulk_modulus_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        pressure: pressure_,
        bulk_modulus: bulk_modulus_,
    })
    return Quantity(result)
