"""
Momentum derivative of kinetic energy is speed
==============================================

The general formula for the kinetic energy of an object features its speed and momentum. This way it can be used
not only in the case of variable mass, but also in the relativistic case.
"""

from sympy import Eq, Derivative
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
)

kinetic_energy = clone_as_function(symbols.kinetic_energy, display_symbol="K(p(v))")
"""
The :symbols:`kinetic_energy` of the object.
"""

momentum = clone_as_function(symbols.momentum, display_symbol="p(v)")
"""
The :symbols:`momentum` of the object.
"""

speed = symbols.speed
"""
The :symbols:`speed` of the object.
"""

law = Eq(
    Derivative(kinetic_energy(momentum(speed)), momentum(speed)),
    speed,
)
"""
:laws:symbol::

:laws:latex::
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
