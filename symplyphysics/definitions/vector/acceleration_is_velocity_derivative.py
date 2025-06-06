"""
Acceleration is velocity derivative
===================================

**Acceleration** is the derivative of velocity w.r.t. time.

**Notes:**

#. Also see :ref:`the scalar counterpart <Acceleration is speed derivative>` of this law.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Acceleration#Instantaneous_acceleration>`__.
"""

from sympy import Eq
from symplyphysics import symbols, Quantity, validate_input, validate_output

from symplyphysics.core.experimental.vectors import clone_as_vector_function, VectorDerivative
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector
from symplyphysics.core.experimental.solvers import solve_for_vector

time = symbols.time
"""
:symbols:`time`.
"""

acceleration = clone_as_vector_function(symbols.acceleration, (time,))
"""
Vector of the body's :symbols:`acceleration` as a function of :attr:`~time`.
"""

velocity = clone_as_vector_function(symbols.speed, (time,))
"""
Vector of the body's velocity as a function of :attr:`~time`. Also see :symbols:`speed`.
"""

law = Eq(acceleration(time), VectorDerivative(velocity(time), time))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    velocity_before_=velocity,
    velocity_after_=velocity,
    time_=time,
)
@validate_output(acceleration)
def calculate_acceleration(
    velocity_before_: QuantityCoordinateVector,
    velocity_after_: QuantityCoordinateVector,
    time_: Quantity,
) -> QuantityCoordinateVector:
    velocity_ = (time / time_) * (velocity_after_ - velocity_before_)

    acceleration_expr = solve_for_vector(law, acceleration(time))
    acceleration_value = acceleration_expr.subs(velocity(time), velocity_).doit()

    return QuantityCoordinateVector.from_expr(acceleration_value)
