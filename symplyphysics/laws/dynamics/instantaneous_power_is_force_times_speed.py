"""
Instantaneous power is force times speed
========================================

*Power* is the rate at which the work is done by a force exerted on an object.
*Instantaneous power* is the power at a specific instant.
"""

from sympy import Eq, cos
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
)
from symplyphysics.core.symbols.quantities import scale_factor

power = symbols.power
"""
The instantaneous :symbols:`power` of :attr:`~force`.
"""

force = symbols.force
"""
The :symbols:`force` exerted on the object.
"""

speed = symbols.speed
"""
The :symbols:`speed` of the object
"""

angle_between_vectors = symbols.angle
r"""
The :symbols:`angle` between the force and velocity vectors.
"""

law = Eq(power, force * speed * cos(angle_between_vectors))
"""
:laws:symbol::

:laws:latex::
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
