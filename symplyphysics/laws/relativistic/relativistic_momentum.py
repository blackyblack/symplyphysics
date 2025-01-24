"""
Relativistic momentum via rest mass and speed
=============================================

Momentum (amount of motion) is a vector physical quantity that is a measure of the mechanical movement
of a body. The relativistic momentum also takes into account speed limits equal to the speed of light.

**Notation:**

#. :quantity_notation:`speed_of_light`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Mass_in_special_relativity#Relativistic_energy%E2%80%93momentum_equation>`__.
"""

from sympy import (Eq, solve, sqrt)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
)

momentum = symbols.momentum
"""
:symbols:`momentum` of the body.
"""

rest_mass = symbols.rest_mass
"""
:symbols:`rest_mass` of the body.
"""

speed = symbols.speed
"""
:symbols:`speed` of the body.
"""

law = Eq(momentum, rest_mass * speed / sqrt(1 - (speed / quantities.speed_of_light)**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(mass_=rest_mass, velocity_=speed)
@validate_output(momentum)
def calculate_momentum(mass_: Quantity, velocity_: Quantity) -> Quantity:
    result_expr = solve(law, momentum, dict=True)[0][momentum]
    result_expr = result_expr.subs({
        rest_mass: mass_,
        speed: velocity_,
    })
    return Quantity(result_expr)
