"""
Momentum is mass times velocity
===============================

Momentum is a physical quantity equal to the product of the object's velocity and its mass.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, symbols)

momentum = Symbol("momentum", units.momentum)
"""
Momentum of the object.

Symbol:
    :code:`p`
"""

velocity = Symbol("velocity", units.velocity)
"""
Velocity of the object.

Symbol:
    :code:`v`
"""

mass = symbols.basic.mass
"""
:attr:`~symplyphysics.symbols.basic.mass` of the object.

Symbol:
    :code:`m`
"""

definition = Eq(momentum, mass * velocity)
"""
:code:`p = m v`

Latex:
    .. math::
        p = m v
"""


@validate_input(velocity_=velocity, mass_=mass)
@validate_output(momentum)
def calculate_momentum(mass_: Quantity, velocity_: Quantity) -> Quantity:
    solved = solve(definition, momentum, dict=True)[0][momentum]
    result_expr = solved.subs({mass: mass_, velocity: velocity_})
    return Quantity(result_expr)
