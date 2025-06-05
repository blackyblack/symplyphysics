"""
Momentum is mass times velocity (Vector)
========================================

An object's *linear momentum* is a vector quantity defined as the product of its mass and velocity vector.

**Links:**

#. `Wikipedia, see first paragraph <https://en.wikipedia.org/wiki/Momentum#>`__.
"""

from sympy import Eq, Expr
from symplyphysics import symbols, Quantity, validate_input, validate_output

from symplyphysics.core.experimental.solvers import solve_for_vector
from symplyphysics.core.experimental.vectors import clone_as_vector_symbol
from symplyphysics.core.experimental.coordinate_systems import (QuantityCoordinateVector,
    combine_coordinate_vectors)

mass = symbols.mass
"""
:symbols:`mass` of the object.
"""

velocity = clone_as_vector_symbol(symbols.speed)
"""
Vector of the object's velocity. See :symbols:`speed`.
"""

momentum = clone_as_vector_symbol(symbols.momentum)
"""
Vector of the object's linear :symbols:`momentum`.
"""

law = Eq(momentum, mass * velocity)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(mass_=mass, velocity_=velocity)
@validate_output(momentum)
def calculate_momentum(
    mass_: Quantity,
    velocity_: QuantityCoordinateVector,
) -> Expr:
    momentum_expr = solve_for_vector(law, momentum)
    momentum_value = momentum_expr.subs({
        mass: mass_,
        velocity: velocity_,
    })

    return combine_coordinate_vectors(momentum_value)


@validate_input(mass_=mass, momentum_=momentum)
@validate_output(velocity)
def calculate_velocity(
    mass_: Quantity,
    momentum_: QuantityCoordinateVector,
) -> Expr:
    velocity_expr = solve_for_vector(law, velocity)
    velocity_value = velocity_expr.subs({
        mass: mass_,
        momentum: momentum_,
    })

    return combine_coordinate_vectors(velocity_value)
