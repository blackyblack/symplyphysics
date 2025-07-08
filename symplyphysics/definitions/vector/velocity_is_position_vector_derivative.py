"""
Velocity is position vector derivative
======================================

Instantaneous **velocity** is the derivative of the body's position vector w.r.t. time.

**Notes:**

#. Also see :ref:`the scalar counterpart <Speed is distance derivative>` of this law.

**Links:**

#. `Wikipedia — Velocity <https://en.wikipedia.org/wiki/Velocity#Instantaneous_velocity>`__.
"""

from sympy import Eq
from symplyphysics import symbols, Quantity, validate_input, validate_output

from symplyphysics.core.experimental.vectors import clone_as_vector_function, VectorDerivative
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

time = symbols.time
"""
:symbols:`time`.
"""

velocity = clone_as_vector_function(symbols.speed, (time,))
"""
Vector of the body's velocity as a function of :attr:`~time`. Also see :symbols:`speed`.
"""

position_vector = clone_as_vector_function(symbols.euclidean_distance, (time,))
"""
Vector of the body's position as a function of :attr:`~time`. Also see :symbols:`euclidean_distance`.
"""

law = Eq(velocity(time), VectorDerivative(position_vector(time), time))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    displacement_=position_vector,
    time_change_=time,
)
@validate_output(velocity)
def calculate_velocity(
    displacement_: QuantityCoordinateVector,
    time_change_: Quantity,
) -> QuantityCoordinateVector:
    # Assuming the velocity is linearly dependent on time.
    position_vector_ = displacement_ * time / time_change_

    result = law.rhs.subs(position_vector(time), position_vector_).doit()

    return QuantityCoordinateVector.from_expr(result)
