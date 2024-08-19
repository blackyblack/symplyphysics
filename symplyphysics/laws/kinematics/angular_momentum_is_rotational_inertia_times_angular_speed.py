"""
Angular momentum is rotational inertia times angular speed
==========================================================

For a rigid body rotating around a fixed axis, the component of its angular momentum parallel
to the rotational axis is found as the product of the body's rotational inertia and the magnitude
of its angular velocity.

**Conditions:**

#. The body is rigid.
#. The axis of rotation is fixed.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    angle_type,
)

angular_momentum = Symbol("angular_momentum", units.length * units.momentum)
"""
Component of the vector of angular momentum parallel to the rotational axis.

Symbol:
    :code:`I`
"""

rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
"""
Rotational inertia of the body.

Symbol:
    :code:`I`
"""

angular_speed = Symbol("angular_speed", angle_type / units.time)
r"""
Angular speed of the body.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

law = Eq(angular_momentum, rotational_inertia * angular_speed)
r"""
:code:`L = I * w`

Latex:
    .. math::
        L = I \omega
"""


@validate_input(rotational_inertia_=rotational_inertia, angular_velocity_=angular_speed)
@validate_output(angular_momentum)
def calculate_angular_momentum(rotational_inertia_: Quantity,
    angular_velocity_: Quantity) -> Quantity:
    result = law.rhs.subs({
        rotational_inertia: rotational_inertia_,
        angular_speed: angular_velocity_,
    })
    return Quantity(result)
