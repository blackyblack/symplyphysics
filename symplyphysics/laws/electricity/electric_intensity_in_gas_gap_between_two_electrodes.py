"""
Electric field in gas gap between two electrodes
================================================

Consider a two-electrode gas-filled gap with flat electrodes. Let one of the electrodes
emit charged particles into a gaseous medium. If the potential collecting these charged
particles is set to another electrode, then a current will flow through the gap.

**Conditions:**

#. Emittivity of the emitter is unlimited, i.e. emission of charged particles is much
   greater than current between electrodes.
#. There is gaseous medium between electrodes.

..
    TODO: find link
"""

from sympy import Eq, solve, sqrt, Rational
from symplyphysics import Quantity, validate_input, validate_output, symbols

electric_field_strength = symbols.electric_field_strength
"""
:symbols:`electric_field_strength` at a point between the electrodes.
"""

coordinate = symbols.position
"""
:symbols:`position` of said point on the axis which runs from the cathode to the anode.
"""

distance = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` between the electrodes.
"""

voltage = symbols.voltage
"""
:symbols:`voltage` between the electrodes.
"""

law = Eq(
    electric_field_strength,
    Rational(3, 2) * (sqrt(coordinate / distance) * (voltage / distance)),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(coordinate_=coordinate,
    distance_between_electrodes_=distance,
    voltage_between_electrodes_=voltage)
@validate_output(electric_field_strength)
def calculate_electric_intensity(coordinate_: Quantity, distance_between_electrodes_: Quantity,
    voltage_between_electrodes_: Quantity) -> Quantity:
    if coordinate_.scale_factor > distance_between_electrodes_.scale_factor:
        raise ValueError("The point of interest should be between electrodes")
    result_expr = solve(law, electric_field_strength, dict=True)[0][electric_field_strength]
    result_expr = result_expr.subs({
        coordinate: coordinate_,
        distance: distance_between_electrodes_,
        voltage: voltage_between_electrodes_,
    })
    return Quantity(result_expr)
