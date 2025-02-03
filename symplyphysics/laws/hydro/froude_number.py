"""
Froude number formula
=====================

The Froude number is based on the speed-to-length ratio as defined by Froude. It has
some analogy with the :ref:`Mach number <Mach number is flow speed over speed of sound>`,
but it not frequently used in the field of theoretical fluid dynamics.

**Notation:**

#. :quantity_notation:`acceleration_due_to_gravity`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Froude_number>`__.

..
    TODO: rename file
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    quantities,
)

flow_speed = symbols.flow_speed
"""
:symbols:`flow_speed`.
"""

characteristic_length = symbols.characteristic_length
"""
:symbols:`characteristic_length` of the fluid container.
"""

froude_number = symbols.froude_number
"""
:symbols:`froude_number`.
"""

law = Eq(froude_number, flow_speed / sqrt(quantities.acceleration_due_to_gravity * characteristic_length))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(velocity_=flow_speed, characteristic_length_=characteristic_length)
@validate_output(froude_number)
def calculate_froude_number(velocity_: Quantity, characteristic_length_: Quantity) -> float:
    result_expr = solve(law, froude_number, dict=True)[0][froude_number]
    result_applied = result_expr.subs({
        flow_speed: velocity_,
        characteristic_length: characteristic_length_
    })
    return convert_to_float(result_applied)
