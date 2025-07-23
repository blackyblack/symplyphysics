"""
Normal force via pressure and vector area
=========================================

Pressure is the proportionality constant between the vector area and the normal force
acting on it.

**Conditions:**

#. Pressure is constant thoughout the given area, which can be achieved by making it
   infinitesimally small.

**Links:**

#. `Wikipedia — Pressure <https://en.wikipedia.org/wiki/Pressure#Formula>`__.

#. `Wikipedia — Normal force <https://en.wikipedia.org/wiki/Normal_force>`__.
"""

from sympy import Eq

from symplyphysics import symbols, Quantity, validate_input, validate_output

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

normal_force = clone_as_vector_symbol(symbols.force, subscript="n")
"""
**Normal force** is the component of a contact force that is normal to the contact surface of the
object.
"""

pressure = symbols.pressure
"""
:symbols:`pressure`.
"""

vector_area = clone_as_vector_symbol(symbols.area)
"""
:ref:`Vector area <Vector area is unit normal times scalar area>`.
"""

law = Eq(normal_force, pressure * vector_area)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    pressure_=pressure,
    vector_area_=vector_area,
)
@validate_output(normal_force)
def calculate_normal_force(
    pressure_: Quantity,
    vector_area_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = law.rhs.subs({
        pressure: pressure_,
        vector_area: vector_area_,
    })

    return QuantityCoordinateVector.from_expr(result)
