"""
Acceleration from force and mass (vector)
=========================================

Newton's second law of motion states that in an inertial frame of reference, the acceleration of
a body in motion is proportional to the net force exerted on it, and its mass is the constant of
proportionality.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Newton's_laws_of_motion#Second_law>`__.
"""

from sympy import Eq, Expr
from symplyphysics import Quantity, validate_input, validate_output, symbols, clone_as_symbol

from symplyphysics.core.experimental.solvers import solve_for_vector
from symplyphysics.core.experimental.vectors import clone_as_vector_symbol
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

mass = clone_as_symbol(symbols.mass, positive=True)
"""
:symbols:`mass` of the object.
"""

acceleration = clone_as_vector_symbol(symbols.acceleration)
"""
:symbols:`acceleration` of the object.
"""

force = clone_as_vector_symbol(symbols.force)
"""
:symbols:`force` exerted on the object.
"""

law = Eq(acceleration, force / mass)
"""
:laws:symbol::

:laws:latex::
"""

# TODO: Derive this law from law of force and momentum
# Condition: mass is constant


@validate_input(mass_=mass, acceleration_=acceleration)
@validate_output(force)
def calculate_force(
    mass_: Quantity,
    acceleration_: QuantityCoordinateVector,
) -> Expr:
    force_expr = solve_for_vector(law, force)
    force_value = force_expr.subs({
        mass: mass_,
        acceleration: acceleration_,
    })

    return QuantityCoordinateVector.from_expr(force_value)


@validate_input(mass_=mass, force_=force)
@validate_output(acceleration)
def calculate_acceleration(
    mass_: Quantity,
    force_: QuantityCoordinateVector,
) -> Expr:
    acceleration_expr = solve_for_vector(law, acceleration)
    acceleration_value = acceleration_expr.subs({
        mass: mass_,
        force: force_,
    })

    return QuantityCoordinateVector.from_expr(acceleration_value)
