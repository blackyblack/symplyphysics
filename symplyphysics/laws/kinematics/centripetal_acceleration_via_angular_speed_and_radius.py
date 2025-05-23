"""
Centripetal acceleration via angular speed and radius
=====================================================

*Centripetal acceleration* is defined as the change in velocity tangential to the velocity vector.

**Links:**

#. `Physics LibreTexts, second part of equation 6.2.5 <https://phys.libretexts.org/Bookshelves/College_Physics/College_Physics_1e_(OpenStax)/06%3A_Uniform_Circular_Motion_and_Gravitation/6.02%3A_Centripetal_Acceleration>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematics import (
    centripetal_acceleration_via_linear_speed_and_radius as centripetal_law,
    speed_via_angular_speed_and_radius as velocities_law,
)

centripetal_acceleration = clone_as_symbol(symbols.acceleration, subscript="n")
"""
Centripetal, or normal, :symbols:`acceleration`.
"""

angular_speed = symbols.angular_speed
"""
:symbols:`angular_speed`.
"""

radius_of_curvature = symbols.radius_of_curvature
"""
Instantaneous :symbols:`radius_of_curvature`.
"""

law = Eq(centripetal_acceleration, angular_speed**2 * radius_of_curvature)
"""
:laws:symbol::

:laws:latex::
"""

# Derive law from expression for linear velocity in circular motion

_centripetal_acceleration_derived = centripetal_law.law.rhs.subs(
    centripetal_law.radius_of_curvature, radius_of_curvature)

_velocities_law_sub = velocities_law.law.subs({
    velocities_law.speed: centripetal_law.speed,
    velocities_law.angular_speed: angular_speed,
    velocities_law.radius_of_curvature: radius_of_curvature,
})

_centripetal_acceleration_derived = solve([
    Eq(centripetal_acceleration, _centripetal_acceleration_derived),
    _velocities_law_sub,
], (centripetal_acceleration, centripetal_law.speed),
    dict=True)[0][centripetal_acceleration]

assert expr_equals(law.rhs, _centripetal_acceleration_derived)


@validate_input(angular_velocity_=angular_speed, curve_radius_=radius_of_curvature)
@validate_output(centripetal_acceleration)
def calculate_centripetal_acceleration(angular_velocity_: Quantity,
    curve_radius_: Quantity) -> Quantity:
    result = law.rhs.subs({
        angular_speed: angular_velocity_,
        radius_of_curvature: curve_radius_,
    })
    return Quantity(result)
