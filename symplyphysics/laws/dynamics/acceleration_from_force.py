"""
Acceleration from force
=======================

Newton's second law of motion states that the acceleration of a body is directly proportional
to the net force exerted on the body.
"""

from sympy import (Eq, solve, sympify)
from symplyphysics import (Vector, Quantity, validate_input, validate_output,
    symbols)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics.vector import acceleration_from_force as acceleration_law_vector

acceleration = symbols.kinematic.acceleration
"""
The :attr:`~symplyphysics.symbols.kinematic.acceleration` of the body.

Symbol:
    a
"""

force = symbols.dynamics.force
"""
The net :attr:`~symplyphysics.symbols.dynamics.force` exerted on the body.

Symbol:
    F
"""

mass = symbols.basic.mass
"""
The :attr:`~symplyphysics.symbols.basic.mass` of the body.

Symbol:
    m
"""

law = Eq(acceleration, force / mass)
r"""
a = F / m

Latex:
    .. math::
        a = \frac{F}{m}
"""

# Derive the same law from vector form

# Scalar law is equivalent to using one-dimensional vectors
_force_vector = Vector([symbols.dynamics.force])
_acceleration_vector = acceleration_law_vector.acceleration_law(_force_vector)
assert len(_acceleration_vector.components) == 1
_acceleration_with_mass = sympify(_acceleration_vector.components[0]).subs(
    acceleration_law_vector.symbols.basic.mass, mass)
assert expr_equals(_acceleration_with_mass, law.rhs)


@validate_input(mass_=mass, acceleration_=acceleration)
@validate_output(force)
def calculate_force(mass_: Quantity, acceleration_: Quantity) -> Quantity:
    result_force_expr = solve(law, force, dict=True)[0][force]
    result_expr = result_force_expr.subs({mass: mass_, acceleration: acceleration_})
    return Quantity(result_expr)
