"""
Voltage is electric field times distance
========================================

**Conditions:**

#. The electric field is constant between the two points. This might be achieved by
   choosing a small enough distance between the points.

**Links:**

#. `Physics LibreTexts, formula 7.3.26 <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_II_-_Thermodynamics_Electricity_and_Magnetism_(OpenStax)/07%3A_Electric_Potential/7.03%3A_Electric_Potential_and_Potential_Difference>`__.
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
