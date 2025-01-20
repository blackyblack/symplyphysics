r"""
Maximum height of body thrown at angle to horizon
=================================================

The maximum height of a projectile above its launch position is a function of
the vertical component of the initial velocity, the angle between the initial
velocity and the horizon, and the magnitude of the free fall acceleration.

**Notation:**

#. :quantity_notation:`acceleration_due_to_gravity`.

**Conditions:**

#. :math:`\mathbf{v}_0 \cdot \mathbf{g} < 0`, i.e. the vector of initial velocity
   :math:`\mathbf{v}_0` must have a non-zero component that is antiparallel to the
   vector of free fall acceleration :math:`\mathbf{g}`.

#. The height is measured with respect to the horizontal plane where the projectile
   was located at initial time.

**Links:**

#. `Physics LibreTexts. Projectile Motion, Maximum Height (3.3.14) <https://phys.libretexts.org/Bookshelves/University_Physics/Physics_(Boundless)/3%3A_Two-Dimensional_Kinematics/3.3%3A_Projectile_Motion>`__.
"""

from sympy import Eq, solve, sin, pi
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
    convert_to_si,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematics import position_via_constant_acceleration_and_time as distance_law
from symplyphysics.laws.kinematics import speed_via_constant_acceleration_and_time as velocity_law
from symplyphysics.laws.geometry import planar_projection_is_cosine as projection_law

height = symbols.height
"""
Maximum :symbols:`height` of the projectile.
"""

initial_speed = clone_as_symbol(symbols.speed, subscript="0")
"""
Initial :symbols:`speed` of the body.
"""

angle = symbols.angle
"""
:symbols:`angle` between the initial velocity and the horizon.
"""

law = Eq(height, initial_speed**2 * sin(angle)**2 / (2 * quantities.acceleration_due_to_gravity))
"""
:laws:symbol::

:laws:latex::
"""

# This law might be derived via "constant_acceleration_movement_is_parabolic" law, "planar_projection_is_cosine" law
# and "accelerated_velocity_from_time" law.

# The law seeks a projection on the horizontal axis, but a projection on the vertical axis is necessary,
# so the angle is represented as a "pi/2 - angle".
_projection_law_applied = projection_law.law.subs({
    projection_law.vector_length: initial_speed,
    projection_law.vector_angle: (pi / 2) - angle,
})
_vertical_projection_derived = solve(_projection_law_applied, projection_law.projection,
    dict=True)[0][projection_law.projection]

# Vertical velocity is zero in the highest point of trajectory.
_velocity_law_applied = velocity_law.law.subs({
    velocity_law.initial_speed: _vertical_projection_derived,
    velocity_law.final_speed: 0,
    velocity_law.acceleration: -quantities.acceleration_due_to_gravity,
})
_time_derived = solve(_velocity_law_applied, velocity_law.time, dict=True)[0][velocity_law.time]

# The acceleration of gravity is directed opposite to the vertical coordinate axis,
## so there is a minus sign before the acceleration.
_height_law_applied = distance_law.law.subs({
    distance_law.initial_speed: _vertical_projection_derived,
    distance_law.time: _time_derived,
    distance_law.acceleration: -quantities.acceleration_due_to_gravity,
    distance_law.initial_position: 0,
})

# Check if derived height is same as declared.
assert expr_equals(_height_law_applied.rhs, law.rhs)


@validate_input(initial_velocity_=initial_speed, angle_=angle)
@validate_output(height)
def calculate_height(initial_velocity_: Quantity, angle_: float | Quantity) -> Quantity:
    result_expr = solve(law, height, dict=True)[0][height]
    result_expr = result_expr.subs({
        initial_speed: initial_velocity_,
        angle: convert_to_si(angle_),
    })
    return Quantity(result_expr)
