"""
Maximum height from initial speed
=================================

The maximum height that a body thrown vertically will rise to depends on the initial speed.

**Links:**

#. `Wikipedia, fifth formula <https://en.wikipedia.org/wiki/Projectile_motion#Maximum_height_of_projectile>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
)

maximum_height = symbols.height
"""
The maximum :symbols:`height` that the object will reach.
"""

initial_speed = symbols.speed
"""
The initial :symbols:`speed` of the object.
"""

law = Eq(maximum_height, initial_speed**2 / (2 * quantities.acceleration_due_to_gravity))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(initial_velocity_=initial_speed)
@validate_output(maximum_height)
def calculate_maximum_height(initial_velocity_: Quantity) -> Quantity:
    result_maximum_height = solve(law, maximum_height, dict=True)[0][maximum_height]
    result_expr = result_maximum_height.subs({initial_speed: initial_velocity_})
    return Quantity(result_expr)
