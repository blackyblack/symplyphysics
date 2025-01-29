"""
Momentum is constant
====================

If there is no external force applied to system of objects, the summary momentum of this
system remains constant during and after any interactions between objects. Summary
momentum of the system is the sum of momentums of every object in this system. Also
applicable for reactive engine simulation.

**Conditions:**

#. The system is closed.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Momentum#Conservation>`__.

..
    TODO: rename file
"""

from sympy import (Derivative, Eq, dsolve)
from symplyphysics import (Quantity, validate_input, validate_output, symbols, clone_as_function)

time = symbols.time
"""
:symbols:`time`.
"""

momentum = clone_as_function(symbols.momentum, [time])
"""
:symbols:`momentum` of the system as a function of :attr:`~time`.
"""

law = Eq(Derivative(momentum(time), time), 0)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(momentum_before_=momentum)
@validate_output(momentum)
def calculate_momentum_after(momentum_before_: Quantity) -> Quantity:
    solved = dsolve(law, momentum(time))
    result_expr = solved.subs("C1", momentum_before_).rhs
    return Quantity(result_expr)
