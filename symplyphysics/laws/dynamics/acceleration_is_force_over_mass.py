"""
Acceleration is force over mass
===============================

Newton's second law of motion states that the acceleration of a body is directly proportional
to the net force exerted on the body.

**Links:**

#. `Britannica <https://www.britannica.com/science/Newtons-laws-of-motion/Newtons-second-law-F-ma>`__.

#. `Wikipedia <https://en.wikipedia.org/wiki/Newton's_laws_of_motion#Second_law>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.laws.dynamics.vector import acceleration_from_force_vector as acceleration_law

from symplyphysics.core.experimental.solvers import solve_for_vector
from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, CoordinateVector

acceleration = symbols.acceleration
"""
The :symbols:`acceleration` of the body.
"""

force = symbols.force
"""
The net :symbols:`force` exerted on the body.
"""

mass = symbols.mass
"""
The :symbols:`mass` of the body.
"""

law = Eq(acceleration, force / mass)
"""
:laws:symbol::

:laws:latex::
"""

# Derive the same law from vector form

# Scalar law is equivalent to using one-dimensional vectors
_force_vector = CoordinateVector([force, 0, 0], CARTESIAN)

_acceleration_vector_derived = solve_for_vector(
    acceleration_law.law,
    acceleration_law.acceleration,
).subs({
    acceleration_law.force: _force_vector,
    acceleration_law.mass: mass,
})

_acceleration_vector_expected = CoordinateVector([law.rhs, 0, 0], CARTESIAN)

assert CoordinateVector.from_expr(_acceleration_vector_derived - _acceleration_vector_expected) == 0


@validate_input(mass_=mass, acceleration_=acceleration)
@validate_output(force)
def calculate_force(mass_: Quantity, acceleration_: Quantity) -> Quantity:
    result_force_expr = solve(law, force, dict=True)[0][force]
    result_expr = result_force_expr.subs({mass: mass_, acceleration: acceleration_})
    return Quantity(result_expr)
