"""
Torque via rotational inertia and angular acceleration
=====================================================

Torque due to a force is a vector physical quantity that describes the effect of the force
on a mechanical object causing its rotational motion.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

torque = symbols.torque
"""
The :symbols:`torque` applied on the object.
"""

rotational_inertia = symbols.rotational_inertia
"""
The :symbols:`rotational_inertia` of the object.
"""

angular_acceleration = symbols.angular_acceleration
r"""
The :symbols:`angular_acceleration` of the object.
"""

law = Eq(torque, rotational_inertia * angular_acceleration)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(moment_of_inertia_=rotational_inertia, angular_acceleration_=angular_acceleration)
@validate_output(torque)
def calculate_moment_of_force(moment_of_inertia_: Quantity,
    angular_acceleration_: Quantity) -> Quantity:
    solved = solve(law, torque, dict=True)[0][torque]
    result_expr = solved.subs({
        rotational_inertia: moment_of_inertia_,
        angular_acceleration: angular_acceleration_
    })
    return Quantity(result_expr)
