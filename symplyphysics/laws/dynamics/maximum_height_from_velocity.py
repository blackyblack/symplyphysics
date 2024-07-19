"""
Maximum height from velocity
============================

The maximum height that a body thrown vertically will rise to depends on the initial velocity.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

maximum_height = Symbol("maximum_height", units.length)
"""
The maximum height that the object will reach.

Symbol:
    h
"""

initial_velocity = Symbol("initial_velocity", units.velocity)
"""
The initial velocity of the object.

Symbol:
    v
"""

law = Eq(maximum_height, initial_velocity**2 / (2 * units.acceleration_due_to_gravity))
r"""
h = v^2 / (2 * g)

Latex:
    :math:`h = \frac{v^2}{2 g}`
"""


@validate_input(initial_velocity_=initial_velocity)
@validate_output(maximum_height)
def calculate_maximum_height(initial_velocity_: Quantity) -> Quantity:
    result_maximum_height = solve(law, maximum_height, dict=True)[0][maximum_height]
    result_expr = result_maximum_height.subs({initial_velocity: initial_velocity_})
    return Quantity(result_expr)
