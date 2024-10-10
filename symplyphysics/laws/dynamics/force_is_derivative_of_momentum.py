"""
Force is derivative of momentum
===============================

Newton's second law of motion can be generalized in terms of linear momentum. Precisely,
the net force exerted on a body is equal to the time derivative of the body's momentum.

**Notes:**

#. Works in relativistic mechanics as well as in classical mechanics.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Momentum#Relation_to_force>`__.
"""

from sympy import Eq, Derivative
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
)

time = symbols.time
"""
:symbols:`time`.
"""

momentum = clone_as_function(symbols.momentum, [time])
"""
The magnitude of the :symbols:`momentum` of the body as a function of :attr:`~time`.
"""

force = clone_as_function(symbols.force, [time])
"""
The magnitude of the net :symbols:`force` exerted on the body as a function of :attr:`~time`.
"""

law = Eq(Derivative(momentum(time), time), force(time))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    momentum_change_=momentum,
    time_change_=time,
)
@validate_output(force)
def calculate_force(
    momentum_change_: Quantity,
    time_change_: Quantity,
) -> Quantity:
    momentum_ = (momentum_change_ / time_change_) * time
    result = law.lhs.subs(momentum(time), momentum_).doit()
    return Quantity(result)
