"""
Mach number is flow speed over speed of sound
=============================================

*Mach number* is a dimensionless quantity in fluid dynamics representing
the ratio of flow speed divided by the speed of sound in the medium.
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
)

mach_number = Symbol("mach_number", dimensionless)
r"""
Local Mach number.

Symbol:
    :code:`M`

Latex:
    :math:`\text{M}`
"""

flow_speed = Symbol("flow_speed", units.velocity)
r"""
Flow speed with respect to internal or external boundaries.

Symbol:
    :code:`u`
"""

speed_of_sound = Symbol("speed_of_sound", units.velocity)
"""
Speed of sound in the medium.

Symbol:
    :code:`c`
"""

law = Eq(mach_number, flow_speed / speed_of_sound)
r"""
:code:`M = u / c`

Latex:
    .. math::
        \text{M} = \frac{u}{c}
"""


@validate_input(velocity_=flow_speed, speed_of_sound_=speed_of_sound)
@validate_output(mach_number)
def calculate_mach_number(velocity_: Quantity, speed_of_sound_: Quantity) -> float:
    result_expr = solve(law, mach_number, dict=True)[0][mach_number]
    result_applied = result_expr.subs({
        flow_speed: velocity_,
        speed_of_sound: speed_of_sound_,
    })
    result = Quantity(result_applied)
    return convert_to_float(result)
