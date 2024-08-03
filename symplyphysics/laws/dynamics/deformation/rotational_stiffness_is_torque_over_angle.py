"""
Rotational stiffness is torque applied over angle
=================================================

*Rotational stiffness* is the extent to which an object resists rotational deformation.
"""

from sympy import Eq
from symplyphysics import (
    units,
    angle_type,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)
from symplyphysics.core.symbols.quantities import scale_factor

rotational_stiffness = Symbol("rotational_stiffness", units.force * units.length / angle_type)
"""
Rotational stiffness of the body.

Symbol:
    :code:`k`
"""

torque = Symbol("torque", units.force * units.length)
r"""
Torque applied to the body.

Symbol:
    :code:`tau`

Latex:
    :math:`\tau`
"""

angular_distance = Symbol("angular_distance", angle_type)
r"""
Angle of deformation.

Symbol:
    :code:`phi`

Latex:
    :math:`\varphi`
"""

law = Eq(rotational_stiffness, torque / angular_distance)
r"""
:code:`k = tau / phi`

Latex:
    .. math::
        k = \frac{\tau}{\varphi}
"""


@validate_input(
    torque_applied_=torque,
    rotation_angle_=angular_distance,
)
@validate_output(rotational_stiffness)
def calculate_rotational_stiffness(
    torque_applied_: Quantity,
    rotation_angle_: float | Quantity,
) -> Quantity:
    result = law.rhs.subs({
        torque: torque_applied_,
        angular_distance: scale_factor(rotation_angle_),
    })
    return Quantity(result)
