"""
Mach number is flow speed over speed of sound
=============================================

*Mach number* is a dimensionless quantity in fluid dynamics representing
the ratio of flow speed divided by the speed of sound in the medium.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Mach_number>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
)

mach_number = symbols.mach_number
"""
Local :symbols:`mach_number`.
"""

flow_speed = symbols.flow_speed
"""
:symbols:`flow_speed` with respect to internal or external boundaries.
"""

speed_of_sound = symbols.speed_of_sound
"""
Local :symbols:`speed_of_sound`.
"""

law = Eq(mach_number, flow_speed / speed_of_sound)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(velocity_=flow_speed, speed_of_sound_=speed_of_sound)
@validate_output(mach_number)
def calculate_mach_number(velocity_: Quantity, speed_of_sound_: Quantity) -> float:
    result_expr = solve(law, mach_number, dict=True)[0][mach_number]
    result_applied = result_expr.subs({
        flow_speed: velocity_,
        speed_of_sound: speed_of_sound_,
    })
    return convert_to_float(result_applied)
