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
    k
"""

torque_applied = Symbol("torque_applied", units.force * units.length)
r"""
Torque applied to the body.

Symbol:
    tau

Latex:
    :math:`\tau`
"""

rotation_angle = Symbol("rotation_angle", angle_type)
r"""
Angle of deformation.

Symbol:
    phi

Latex:
    :math:`\varphi`
"""

law = Eq(rotational_stiffness, torque_applied / rotation_angle)
r"""
k = tau / phi

Latex:
    :math:`k = \frac{\tau}{\varphi}`
"""

@validate_input(
    torque_applied_=torque_applied,
    rotation_angle_=rotation_angle,
)
@validate_output(rotational_stiffness)
def calculate_rotational_stiffness(
    torque_applied_: Quantity,
    rotation_angle_: float | Quantity,
) -> Quantity:
    result = law.rhs.subs({
        torque_applied: torque_applied_,
        rotation_angle: scale_factor(rotation_angle_),
    })
    return Quantity(result)
