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
from symplyphysics import Quantity, validate_input, validate_output, symbols

angular_momentum = symbols.angular_momentum
"""
Component of the vector of :symbols:`angular_momentum` parallel to the rotational axis.
"""

rotational_inertia = symbols.rotational_inertia
"""
:symbols:`rotational_inertia` of the body.
"""

angular_speed = symbols.angular_speed
"""
:symbols:`angular_speed` of the body.
"""

law = Eq(angular_momentum, rotational_inertia * angular_speed)
"""
:laws:symbol::

:laws:latex::
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
