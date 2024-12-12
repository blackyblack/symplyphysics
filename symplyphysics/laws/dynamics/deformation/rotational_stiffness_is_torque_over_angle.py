"""
Rotational stiffness is torque applied over angle
=================================================

*Rotational stiffness* is the extent to which an object resists rotational deformation.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Stiffness#Rotational_stiffness>`__.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.symbols.quantities import scale_factor

rotational_stiffness = symbols.torsion_stiffness
"""
Rotational stiffness of the body. See :symbols:`torsion_stiffness`.
"""

torque = symbols.torque
"""
:symbols:`torque` applied to the body.
"""

angular_distance = symbols.angular_distance
"""
Angle of deformation. See :symbols:`angular_distance`
"""

law = Eq(rotational_stiffness, torque / angular_distance)
"""
:laws:symbol::

:laws:latex::
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
