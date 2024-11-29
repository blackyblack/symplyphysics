"""
Time of flight of a projectile via maximum height
=================================================

Let the body be thrown horizontally with some initial velocity. Then the time of the motion
of this body until it reaches the initial elevation depends only on the height and acceleration
of free fall.

**Conditions:**

#. The acceleration of gravity stays constant throughout the movement of the body.

#. The initial and final height of the body with respect to the "ground" are the same,
   i.e. the body lands on a point with the same vertical coordinate as at the start
   of the flight.

**Links:**

#. `Physics LibreTexts. Zero Launch Angle, Duration of Flight (3.3.22) <https://phys.libretexts.org/Bookshelves/University_Physics/Physics_(Boundless)/3%3A_Two-Dimensional_Kinematics/3.3%3A_Projectile_Motion>`__.

..
    TODO rename file
"""

from sympy import (Eq, solve, sqrt)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    quantities,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematics import position_via_constant_acceleration_and_time as distance_law

time = symbols.time
"""
:symbols:`time` of flight of the projectile.
"""

height = symbols.height
"""
Maximum :symbols:`height` which the projectile reaches during its motion.
"""

law = Eq(time, sqrt(2 * height / quantities.acceleration_due_to_gravity))
"""
:laws:symbol::

:laws:latex::
"""

# This law might be derived via "constant_acceleration_movement_is_parabolic" law.
# Horizontal vector of movement does not change falling time.

_distance_law_applied = distance_law.law.subs({
    distance_law.initial_speed: 0,
    distance_law.acceleration: quantities.acceleration_due_to_gravity,
    distance_law.initial_position: 0,
    distance_law.final_position: height,
})
_time_derived = solve(_distance_law_applied, distance_law.time, dict=True)[1][distance_law.time]

# Check if derived movement time is same as declared.
assert expr_equals(_time_derived, law.rhs)


@validate_input(height_=height)
@validate_output(time)
def calculate_movement_time(height_: Quantity) -> Quantity:
    result_expr = solve(law, time, dict=True)[0][time]
    result_expr = result_expr.subs({
        height: height_,
    })
    return Quantity(result_expr)
