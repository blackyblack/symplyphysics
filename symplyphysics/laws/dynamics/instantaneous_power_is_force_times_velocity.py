"""
Instantaneous power is force times speed
========================================

*Power* is the rate at which the work is done by a force exerted on an object.
*Instantaneous power* is the power at a specific instant.
"""

from sympy import Eq, cos
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    angle_type,
)
from symplyphysics.core.symbols.quantities import scale_factor

power = Symbol("power", units.power)
"""
The instantaneous power of force :math:`F`.

Symbol:
    :code:`P`
"""

force = symbols.dynamics.force
"""
The :attr:`~symplyphysics.symbols.dynamics.force` exerted on the object.

Symbol:
    :code:`F`
"""

speed = Symbol("speed", units.speed)
"""
The speed of the object

Symbol:
    :code:`v`
"""

angle_between_vectors = Symbol("angle_between_vectors", angle_type)
r"""
The angle between the force and velocity vectors.

Symbol:
    :code;`phi`

Latex:
    :math:`\varphi`
"""

law = Eq(power, force * speed * cos(angle_between_vectors))
r"""
:code:`P = F * v * cos(phi)`

Latex:
    .. math::
        P = F v \cos{\varphi}
"""


# TODO: derive law from [its vector counterpart](./vector/instantaneous_power_is_force_dot_velocity.py)


@validate_input(force_=force, speed_=speed, angle_=angle_between_vectors)
@validate_output(power)
def calculate_power(force_: Quantity, speed_: Quantity, angle_: Quantity | float) -> Quantity:
    angle_ = scale_factor(angle_)
    result = law.rhs.subs({
        force: force_,
        speed: speed_,
        angle_between_vectors: angle_,
    })
    return Quantity(result)
