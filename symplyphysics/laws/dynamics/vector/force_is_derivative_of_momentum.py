"""
Force is derivative of momentum (vector)
========================================

Newton's second law of motion can be generalized in terms of linear momentum. Precisely,
the net force exerted on a body is equal to the time derivative of the body's momentum.

**Notes:**

#. Works in relativistic mechanics as well as in classical mechanics.

#. See :ref:`scalar counterpart <Force is derivative of momentum>` of this law.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Momentum#Relation_to_force>`__.
"""

from sympy import Eq, Expr

from symplyphysics import Quantity, symbols, validate_input, validate_output

from symplyphysics.core.experimental.coordinate_systems import (QuantityCoordinateVector,
    combine_coordinate_vectors)
from symplyphysics.core.experimental.vectors import VectorDerivative, clone_as_vector_function
from symplyphysics.core.experimental.solvers import solve_for_vector

time = symbols.time
"""
:symbols:`time`.
"""

momentum = clone_as_vector_function(symbols.momentum, (time,))
"""
The magnitude of the :symbols:`momentum` of the body as a function of :attr:`~time`.
"""

force = clone_as_vector_function(symbols.force, (time,))
"""
The magnitude of the net :symbols:`force` exerted on the body as a function of :attr:`~time`.
"""

force_law = Eq(force(time), VectorDerivative(momentum(time), time))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    momentum_before_=momentum,
    momentum_after_=momentum,
    time_=time,
)
@validate_output(force)
def calculate_force(
    momentum_before_: QuantityCoordinateVector,
    momentum_after_: QuantityCoordinateVector,
    time_: Quantity,
) -> Expr:
    momentum_ = (time / time_) * (momentum_after_ - momentum_before_)

    solved = solve_for_vector(force_law, force(time))
    force_ = solved.subs(momentum(time), momentum_).doit()
    force_ = combine_coordinate_vectors(force_)

    return force_
