"""
Maximum height from initial speed
=================================

The maximum height that a body thrown vertically will rise to depends on the initial speed.

..
    TODO Rename file
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

maximum_height = Symbol("maximum_height", units.length)
"""
The maximum height that the object will reach.

Symbol:
    :code:`h`
"""

initial_speed = Symbol("initial_speed", units.velocity)
"""
The initial speed of the object.

Symbol:
    :code:`v`
"""

law = Eq(maximum_height, initial_speed**2 / (2 * units.acceleration_due_to_gravity))
r"""
:code:`h = v^2 / (2 * g)`

Latex:
    .. math::
        h = \frac{v^2}{2 g}
"""


@validate_input(initial_velocity_=initial_speed)
@validate_output(maximum_height)
def calculate_maximum_height(initial_velocity_: Quantity) -> Quantity:
    result_maximum_height = solve(law, maximum_height, dict=True)[0][maximum_height]
    result_expr = result_maximum_height.subs({initial_speed: initial_velocity_})
    return Quantity(result_expr)
