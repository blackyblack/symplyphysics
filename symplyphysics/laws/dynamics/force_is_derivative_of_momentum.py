"""
Force is derivative of momentum
===============================

Newton's second law of motion can be generalized in terms of linear momentum. Precisely,
the net force exerted on a body is equal to the time derivative of the body's momentum.

**Notes:**

#. Works in relativistic mechanics as well as in classical mechanics.
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

momentum = Function("momentum", units.momentum)
"""
The magnitude of the momentum of the body.

Symbol:
    p
"""

force = Function("force", units.force)
"""
The magnitude of the net force exerted on the body.

Symbol:
    F
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    t
"""

law = Eq(Derivative(momentum(time), time), force(time))
r"""
dp/dt = F(t)

Latex:
    :math:`\frac{d p}{d t} = F(t)`
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
