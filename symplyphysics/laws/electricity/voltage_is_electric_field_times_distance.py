"""
Voltage is electric field times distance
========================================

**Conditions:**

#. The electric field is constant between the two points. This might be achieved by
   choosing a small enough distance between the points.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

voltage = symbols.voltage
"""
:symbols:`voltage` between two points.
"""

electric_field_strength = symbols.electric_field_strength
"""
:symbols:`electric_field_strength`.
"""

distance = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` between two points.
"""

law = Eq(voltage, electric_field_strength * distance)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    electric_field_strength_=electric_field_strength,
    distance_=distance,
)
@validate_output(voltage)
def calculate_voltage(
    electric_field_strength_: Quantity,
    distance_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        electric_field_strength: electric_field_strength_,
        distance: distance_,
    })
    return Quantity(result)
