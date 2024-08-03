"""
Torque via rotational inertia and angular acceleration
=====================================================

Torque due to a force is a vector physical quantity that describes the effect of the force
on a mechanical object causing its rotational motion.

..
    TODO Rename file
"""

from sympy import Eq, solve
from symplyphysics import (angle_type, units, Quantity, Symbol, validate_input, validate_output)

torque = Symbol("torque", units.force * units.length)
r"""
The torque applied on the object.

Symbol:
    :code:`tau`

Latex:
    :math:`\tau`
"""

rotational_inertia = Symbol("rotational_inertia", units.mass * units.area)
"""
The rotational inertia of the object.

Symbol:
    :code:`I`
"""

angular_acceleration = Symbol("angular_acceleration", angle_type / (units.time**2))
r"""
The angular acceleration of the object.

Symbol:
    :code:`epsilon`

Latex:
    :math:`\varepsilon`
"""

law = Eq(torque, rotational_inertia * angular_acceleration)
r"""
:code:`tau = I * epsilon`

Latex:
    .. math::
        \tau = I \varepsilon
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
