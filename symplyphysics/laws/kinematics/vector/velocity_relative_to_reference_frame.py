"""
Velocity relative to reference frame
====================================

For any reference frame, whether it is inertial or not, the motion relative to it can be described
using the position vector relative to that frame's origin.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Velocity#Instantaneous_velocity>`__.

#. `Wikipedia (ru) <https://ru.wikipedia.org/wiki/%D0%A2%D0%B5%D0%BE%D1%80%D0%B5%D0%BC%D0%B0_%D0%BE_%D1%81%D0%BB%D0%BE%D0%B6%D0%B5%D0%BD%D0%B8%D0%B8_%D1%81%D0%BA%D0%BE%D1%80%D0%BE%D1%81%D1%82%D0%B5%D0%B9#%D0%9E%D0%B1%D1%81%D1%83%D0%B6%D0%B4%D0%B5%D0%BD%D0%B8%D0%B5>`__.
"""

from sympy import Eq
from symplyphysics import symbols, validate_input, validate_output, Quantity

from symplyphysics.core.experimental.vectors import clone_as_vector_function, SymDerivative
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

time = symbols.time
"""
:symbols:`time`.
"""

relative_velocity = clone_as_vector_function(symbols.speed, (time,))
"""
Vector of the body's speed relative to :math:`S` as a function of :attr:`~time`.
"""

position_vector = clone_as_vector_function(symbols.distance_to_origin, (time,))
"""
Vector of the body's position relative to :math:`S` as a function of :attr:`~time`.
"""

law = Eq(relative_velocity(time), SymDerivative(position_vector(time), time))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    position_before_=position_vector,
    position_after_=position_vector,
    time_change_=time,
)
@validate_output(relative_velocity)
def calculate_relative_velocity(
    position_before_: QuantityCoordinateVector,
    position_after_: QuantityCoordinateVector,
    time_change_: Quantity,
) -> QuantityCoordinateVector:
    position_vector_ = (time / time_change_) * (position_after_ - position_before_)

    result = law.rhs.subs(position_vector(time), position_vector_).doit()

    return QuantityCoordinateVector.from_expr(result)
