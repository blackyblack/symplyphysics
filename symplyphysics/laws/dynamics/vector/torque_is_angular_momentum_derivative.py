"""
Torque is angular momentum derivative
=====================================

In the case of a single particle, the total vector sum of all the external torques acting on a
particle is equal to the time rate change of the angular momentum of that particle.

In the case of a system of particles, the net external torque acting on a system of particles is
equal to the time rate change of the system's total angular momentum.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Torque#Relationship_with_the_angular_momentum>`__.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols

from symplyphysics.core.experimental.vectors import clone_as_vector_function, VectorDerivative
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector
from symplyphysics.core.experimental.solvers import solve_for_vector

time = symbols.time
"""
:symbols:`time`.
"""

angular_momentum = clone_as_vector_function(symbols.angular_momentum, (time,))
"""
Pseudovector of net :symbols:`angular_momentum` as a function of :attr:`~time`.
"""

torque = clone_as_vector_function(symbols.torque, (time,))
"""
Pseudovector of net :symbols:`torque` as a function of :attr:`~time`.
"""

law = Eq(torque(time), VectorDerivative(angular_momentum(time), time))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    angular_momentum_before_=angular_momentum,
    angular_momentum_after_=angular_momentum,
    time_=time,
)
@validate_output(torque)
def calculate_torque(
    angular_momentum_before_: QuantityCoordinateVector,
    angular_momentum_after_: QuantityCoordinateVector,
    time_: Quantity,
) -> QuantityCoordinateVector:
    angular_momentum_ = (time / time_) * (angular_momentum_after_ - angular_momentum_before_)

    result = solve_for_vector(
        law.subs(angular_momentum(time), angular_momentum_),
        torque(time),
    ).doit()

    return QuantityCoordinateVector.from_expr(result)
