"""
Time of flight of a projectile via initial velocity
===================================================

Let the body be thrown horizontally with some initial velocity. Then the time of the motion
of this body until it reaches the initial elevation depends only initial velocity and acceleration
of free fall.

**Conditions:**

#. The acceleration of gravity stays constant throughout the movement of the body.

#. The initial and final height of the body with respect to the "ground" are the same,
   i.e. the body lands on a point with the same vertical coordinate as at the start
   of the flight.

**Links:**

#. `Physics LibreTexts. Projectile Motion, Maximum Height (3.3.13) <https://phys.libretexts.org/Bookshelves/University_Physics/Physics_(Boundless)/3%3A_Two-Dimensional_Kinematics/3.3%3A_Projectile_Motion>`__.

..
    TODO rename file
"""

from sympy import Eq, solve, sin, pi
from symplyphysics import Quantity, validate_input, validate_output, symbols, quantities
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematics import position_via_constant_acceleration_and_time as distance_law
from symplyphysics.laws.geometry import planar_projection_is_cosine as projection_law

time = symbols.time
"""
:symbols:`time` of flight of the projectile.
"""

initial_speed = symbols.speed
"""
Initial :symbols:`speed` of the projectile.
"""

angle = symbols.angle
"""
:symbols:`angle`
"""

law = Eq(time, 2 * initial_speed * sin(angle) / quantities.acceleration_due_to_gravity)
"""
:laws:symbol::

:laws:latex::
"""

# This law might be derived via "constant_acceleration_movement_is_parabolic" law
# and "planar_projection_is_cosine" law.

# The law seeks a projection on the horizontal axis, but a projection on the vertical axis is necessary,
# so the angle is represented as a "pi/2 - angle".
_projection_law_applied = projection_law.law.subs({
    projection_law.vector_length: initial_speed,
    projection_law.vector_angle: (pi / 2) - angle,
})
_projection_derived = solve(_projection_law_applied, projection_law.projection,
    dict=True)[0][projection_law.projection]

# The acceleration of gravity is directed opposite to the vertical coordinate axis,
## so there is a minus sign before the acceleration.
_distance_law_applied = distance_law.law.subs({
    distance_law.initial_speed: _projection_derived,
    distance_law.acceleration: -quantities.acceleration_due_to_gravity,
    distance_law.initial_position: 0,
    distance_law.final_position: 0,
})
_time_derived = solve(_distance_law_applied, distance_law.time, dict=True)[1][distance_law.time]

# Check if derived movement time is same as declared.
assert expr_equals(_time_derived, law.rhs)


@validate_input(initial_velocity_=initial_speed, angle_=angle)
@validate_output(time)
def calculate_movement_time(initial_velocity_: Quantity, angle_: float | Quantity) -> Quantity:
    result_expr = solve(law, time, dict=True)[0][time]
    result_expr = result_expr.subs({
        initial_speed: initial_velocity_,
        angle: angle_,
    })
    return Quantity(result_expr)
