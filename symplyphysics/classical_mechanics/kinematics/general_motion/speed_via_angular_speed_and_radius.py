"""
Speed via angular speed and radius
==================================

Speed of a rotating body can be calculated using its angular speed and instantaneous
radius of curvature of the body's path.

**Links:**

#. `Physics LibreTexts, first part of formula 6.1.9 <https://phys.libretexts.org/Bookshelves/College_Physics/College_Physics_1e_(OpenStax)/06%3A_Uniform_Circular_Motion_and_Gravitation/6.01%3A_Rotation_Angle_and_Angular_Velocity>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols

speed = symbols.speed
"""
Linear :symbols:`speed`.
"""

angular_speed = symbols.angular_speed
"""
:symbols:`angular_speed`.
"""

radius_of_curvature = symbols.radius_of_curvature
"""
Instantaneous :symbols:`radius_of_curvature`.
"""

law = Eq(speed, angular_speed * radius_of_curvature)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(angular_velocity_=angular_speed, curve_radius_=radius_of_curvature)
@validate_output(speed)
def calculate_linear_velocity(angular_velocity_: Quantity, curve_radius_: Quantity) -> Quantity:
    solved = solve(law, speed, dict=True)[0][speed]
    result_expr = solved.subs({
        angular_speed: angular_velocity_,
        radius_of_curvature: curve_radius_
    })
    return Quantity(result_expr)
