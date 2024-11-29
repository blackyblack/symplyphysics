"""
Horizontal displacement of projectile
=====================================

Let's say we throw the body at an angle to the horizon with some initial velocity.
Then the horizontal displacement, or **range**, of the body depends on the initial
speed, the angle of the throw and the acceleration of free fall.

**Conditions:**

#. The acceleration due to gravity is constant throughout the movement
   of the projectile.

#. The elevation of the initial and final points of the projectile are
   the same.

**Links:**

#. `Physics LibreTexts. Projectile Motion, Range (3.3.15) <https://phys.libretexts.org/Bookshelves/University_Physics/Physics_(Boundless)/3%3A_Two-Dimensional_Kinematics/3.3%3A_Projectile_Motion>`__.
"""

from sympy import Eq, solve, sin, simplify
from symplyphysics import Quantity, validate_input, validate_output, symbols, quantities
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematics import position_via_constant_speed_and_time as distance_law
from symplyphysics.laws.geometry import planar_projection_is_cosine as projection_law
from symplyphysics.laws.gravity import maximum_movement_time_of_a_body_thrown_at_an_angle_to_horizon as time_law

horizontal_displacement = symbols.euclidean_distance
"""
Horizontal displacement, or **range**, of the projectile. See :symbols:`euclidean_distance`.
"""

initial_speed = symbols.speed
"""
Initial :symbols:`speed` of the projectile.
"""

angle = symbols.angle
"""
:symbols:`angle` between the initial velocity and the horizon.
"""

law = Eq(horizontal_displacement, initial_speed**2 * sin(2 * angle) / quantities.acceleration_due_to_gravity)
"""
:laws:symbol::

:laws:latex::
"""

# This law might be derived via "distance_from_constant_velocity" law, "planar_projection_is_cosine" law
# and "maximum_movement_time_of_a_body_thrown_at_an_angle_to_horizon" law.

_projection_law_applied = projection_law.law.subs({
    projection_law.vector_length: initial_speed,
    projection_law.vector_angle: angle,
})
_horizontal_projection_derived = solve(_projection_law_applied, projection_law.projection,
    dict=True)[0][projection_law.projection]

_time_law_applied = time_law.law.subs({
    time_law.initial_speed: initial_speed,
    time_law.angle: angle,
})
_time_derived = solve(_time_law_applied, time_law.movement_time, dict=True)[0][time_law.movement_time]

_range_law_applied = distance_law.law.subs({
    distance_law.initial_position: 0,
    distance_law.speed: _horizontal_projection_derived,
    distance_law.time: _time_derived,
})

# Check if derived range is same as declared.
assert expr_equals(simplify(_range_law_applied.rhs), law.rhs)


@validate_input(initial_velocity_=initial_speed, angle_=angle)
@validate_output(horizontal_displacement)
def calculate_range(
    initial_velocity_: Quantity,
    angle_: float | Quantity,
) -> Quantity:
    result_expr = solve(law, horizontal_displacement, dict=True)[0][horizontal_displacement]
    result_expr = result_expr.subs({
        initial_speed: initial_velocity_,
        angle: angle_,
    })
    return Quantity(result_expr)
