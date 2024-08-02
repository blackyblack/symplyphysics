"""
Elastic energy density of bulk compression via pressure
=======================================================

Volumetric density of the elastic energy of a body that is under bulk
compression is proportional to the square of the pressure in the body
and inversely proportional to the bulk modulus of the body's material.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

elastic_energy_density = Symbol("elastic_energy_density", units.energy / units.volume)
"""
Elastic energy of the deformed body per unit of its volume.

Symbol:
    :code:`u`
"""

pressure = Symbol("pressure", units.pressure)
"""
Pressure in the deformed body.

Symbol:
    :code:`P`
"""

bulk_modulus = Symbol("bulk_modulus", units.pressure)
"""
Bulk modulus of the material.

Symbol:
    :code:`K`
"""

law = Eq(elastic_energy_density, pressure**2 / (2 * bulk_modulus))
r"""
:code:`u = P^2 / (2 * K)`

Latex:
    .. math::
        u = \frac{P^2}{2 K}
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
