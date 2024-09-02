"""
Speed via angular speed and radius
==================================

Speed of a rotating body can be calculated using its angular speed and instantaneous
radius of curvature of the body's path.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, angle_type, validate_input, validate_output)

speed = Symbol("speed", units.velocity)
"""
Linear speed.

Symbol:
    :code:`v`
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
r"""
Instantaneous radius of curvature.

Symbol:
    :code:`r`
"""

law = Eq(speed, angular_speed * radius_of_curvature)
r"""
:code:`v = w * r`

Latex:
    .. math::
        v = \omega r
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
