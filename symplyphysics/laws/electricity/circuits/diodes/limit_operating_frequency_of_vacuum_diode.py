"""
Limit operating frequency of vacuum diode
=========================================

The vacuum diode has a limit operating frequency. At frequencies higher, the diode stops
functioning normally, since the electrons will not have time to move from the cathode to
the anode.

**Notation:**

#. :quantity_notation:`elementary_charge`.
#. :quantity_notation:`electron_rest_mass`.

..
    TODO: find link
"""

from sympy import Eq, solve, sqrt
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.quantities import elementary_charge, electron_rest_mass

limit_frequency = symbols.temporal_frequency
"""
Limit operating :symbols:`temporal_frequency` of the vacuum diode.
"""

distance = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` between electrodes.
"""

voltage = symbols.voltage
"""
:symbols:`voltage` between cathode and anode.
"""

law = Eq(
    limit_frequency,
    sqrt(2 * elementary_charge * voltage / electron_rest_mass) / (6 * distance))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(distance_between_electrodes_=distance, voltage_=voltage)
@validate_output(limit_frequency)
def calculate_limit_operating_frequency(distance_between_electrodes_: Quantity,
    voltage_: Quantity) -> Quantity:
    result_expr = solve(law, limit_frequency, dict=True)[0][limit_frequency]
    result_expr = result_expr.subs({
        distance: distance_between_electrodes_,
        voltage: voltage_,
    })
    return Quantity(result_expr)
