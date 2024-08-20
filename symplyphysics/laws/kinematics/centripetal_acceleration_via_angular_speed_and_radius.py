"""
Centripetal acceleration via angular speed and radius
=====================================================

*Centripetal acceleration* is defined as the change in velocity tangential to the velocity vector.
"""

from sympy import Eq, solve
from symplyphysics import (
    clone_symbol,
    symbols,
    units,
    Quantity,
    Symbol,
    angle_type,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematics import (
    centripetal_acceleration_via_linear_speed_and_radius as centripetal_law,
    speed_via_angular_speed_and_radius as velocities_law,
)

centripetal_acceleration = clone_symbol(symbols.kinematics.acceleration, "centripetal_acceleration")
r"""
Centripetal, or normal, :attr:`~symplyphysics.symbols.kinematics.acceleration`.

Symbol:
    :code:`a_n`

Latex:
    :math:`a_n`
"""

angular_speed = Symbol("angular_speed", angle_type / units.time)
r"""
Angular speed.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

radius_of_curvature = Symbol("radius_of_curvature", units.length)
"""
Instantaneous radius of curvature.

Symbol:
    :code:`r`
"""

law = Eq(centripetal_acceleration, angular_speed**2 * radius_of_curvature)
r"""
:code:`a_n = w^2 * r`

Latex:
    .. math::
        a_n = \omega^2 r
"""

# Derive law from expression for linear velocity in circular motion

_centripetal_acceleration_derived = centripetal_law.law.rhs.subs(centripetal_law.radius_of_curvature,
    radius_of_curvature)

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
