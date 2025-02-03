"""
Nusselt number formula
======================

Nusselt number is the ratio of total heat transfer to conductive heat transfer at a
boundary in a fluid. It can be expressed using the heat transfer coefficient of the
flow, characteristic length of the system, and thermal conductivity of the fluid.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Nusselt_number#Definition>`__.

..
    TODO: rename file
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols
)

heat_transfer_coefficient = symbols.heat_transfer_coefficient
"""
:symbols:`heat_transfer_coefficient` of the flow.
"""

characteristic_length = symbols.characteristic_length
"""
:symbols:`characteristic_length` of the system.
"""

thermal_conductivity = symbols.thermal_conductivity
"""
:symbols:`thermal_conductivity` of the fluid.
"""

nusselt_number = symbols.nusselt_number
"""
:symbols:`nusselt_number`.
"""

law = Eq(nusselt_number, heat_transfer_coefficient * characteristic_length / thermal_conductivity)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    heat_transfer_coefficient_=heat_transfer_coefficient,
    characteristic_length_=characteristic_length,
    thermal_conductivity_=thermal_conductivity,
)
@validate_output(nusselt_number)
def calculate_nusselt_number(heat_transfer_coefficient_: Quantity, characteristic_length_: Quantity,
    thermal_conductivity_: Quantity) -> float:
    result_expr = solve(law, nusselt_number, dict=True)[0][nusselt_number]
    result_applied = result_expr.subs({
        heat_transfer_coefficient: heat_transfer_coefficient_,
        characteristic_length: characteristic_length_,
        thermal_conductivity: thermal_conductivity_
    })
    return convert_to_float(result_applied)
