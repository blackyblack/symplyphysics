"""
Momentum derivative of kinetic energy is speed
==============================================

The general formula for the kinetic energy of an object features its speed and momentum. This way it can be used
not only in the case of variable mass, but also in the relativistic case.
"""

from sympy import Eq, Derivative
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
)

kinetic_energy = Function("kinetic_energy", units.energy)
r"""
The kinetic energy of the object.

Symbol:
    :code:`K(p(v))`
"""

momentum = Function("momentum", units.momentum)
"""
The momentum of the object.

Symbol:
    :code:`p(v)`
"""

speed = Symbol("speed", units.velocity)
"""
The speed of the object.

Symbol:
    :code:`v`
"""

law = Eq(
    Derivative(kinetic_energy(momentum(speed)), momentum(speed)),
    speed,
)
r"""
:code:`Derivative(K(p(v), p(v))) = v`

Latex:
    .. math::
        \frac{d K(p(v))}{d p(v)} = v
"""

# TODO: derive from the differential definition of work and the generalized Newton's second law


@validate_input(
    kinetic_energy_change_=kinetic_energy,
    momentum_change_=momentum,
)
@validate_output(speed)
def calculate_speed(
    kinetic_energy_change_: Quantity,
    momentum_change_: Quantity,
) -> Quantity:
    kinetic_energy_ = (kinetic_energy_change_ / momentum_change_) * momentum(speed)
    result = law.lhs.subs(
        kinetic_energy(momentum(speed)),
        kinetic_energy_,
    ).doit()
    return Quantity(result)
