"""
Damping force is proportional to velocity
=========================================

Damping force is an external (relative to an object) force that drains energy from the object,
reducing the motion of the object. It is a model used, for example, to describe the motion
of an oscillator.

**Links:**

#. `Physics LibreTexts, formula 8.3.1 <https://phys.libretexts.org/Courses/University_of_California_Davis/UCD%3A_Physics_9HA__Classical_Mechanics/8%3A_Small_Oscillations/8.3%3A_Damping_and_Resonance>`__.
"""

from sympy import Eq
from symplyphysics import symbols, Quantity, validate_input, validate_output

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector
from symplyphysics.core.experimental.solvers import solve_for_vector

damping_constant = symbols.damping_constant
"""
Non-negative :symbols:`damping_constant`.
"""

velocity = clone_as_vector_symbol(symbols.speed)
"""
Vector of the body's velocity. See :symbols:`speed`.
"""

damping_force = clone_as_vector_symbol(symbols.force)
"""
Vector of the damping :symbols:`force` exerted on the body.
"""

law = Eq(damping_force, -1 * damping_constant * velocity)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(damping_constant_=damping_constant, velocity_=velocity)
@validate_output(damping_force)
def calculate_damping_force(
    damping_constant_: Quantity,
    velocity_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    force_expr = solve_for_vector(law, damping_force)
    force_value = force_expr.subs({
        damping_constant: damping_constant_,
        velocity: velocity_,
    })

    return QuantityCoordinateVector.from_expr(force_value)


@validate_input(damping_constant_=damping_constant, damping_force_=damping_force)
@validate_output(velocity)
def calculate_velocity(
    damping_constant_: Quantity,
    damping_force_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    velocity_expr = solve_for_vector(law, velocity)
    velocity_value = velocity_expr.subs({
        damping_constant: damping_constant_,
        damping_force: damping_force_,
    })

    return QuantityCoordinateVector.from_expr(velocity_value)


# UNIQUE_LAW_ID: 783
