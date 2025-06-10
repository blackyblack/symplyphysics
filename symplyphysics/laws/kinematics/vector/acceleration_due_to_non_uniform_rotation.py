"""
Acceleration due to non-uniform rotation
========================================

Imagine two reference frames, one of which is fixed (:math:`S`) and the other is moving
(:math:`S'`). When :math:`S'` rotates around :math:`S` in a non-uniform way, the acceleration of
some body :math:`B` in :math:`S` has a component corresponding to that non-uniform rotation of
:math:`S'`. It is part of the transfer acceleration of body :math:`B` in :math:`S`.

**Notation:**

#. :math:`\\left[ \\vec a, \\vec b \\right]` (:code:`cross(a, b)`) is vector product of
   :math:`\\vec a` and :math:`\\vec b`.

**Links:**

#. `Wikipedia <https://ru.wikipedia.org/wiki/%D0%A1%D0%BB%D0%BE%D0%B6%D0%BD%D0%BE%D0%B5_%D0%B4%D0%B2%D0%B8%D0%B6%D0%B5%D0%BD%D0%B8%D0%B5#%D0%A3%D1%81%D0%BA%D0%BE%D1%80%D0%B5%D0%BD%D0%B8%D0%B5>`__.

..
    TODO find English link
"""

from sympy import Eq
from symplyphysics import validate_input, validate_output, Quantity, symbols

from symplyphysics.core.experimental.vectors import (clone_as_vector_symbol,
    clone_as_vector_function, VectorCross, VectorDerivative)
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

non_uniform_rotation_acceleration = clone_as_vector_symbol(
    symbols.acceleration,
    display_symbol="a_rot",
    display_latex="{\\vec a}_\\text{rot}",
)
"""
TODO
"""

time = symbols.time
"""
:symbols:`time`.
"""

angular_velocity = clone_as_vector_function(symbols.angular_speed, (time,))
"""
Pseudovector of the body's angular velocity as a function of :attr:`~time`. See
:symbols:`angular_speed`.
"""

position_vector = clone_as_vector_symbol(symbols.distance_to_origin)
"""
The body's position vector. See :symbols:`distance_to_origin`.
"""

law = Eq(
    non_uniform_rotation_acceleration,
    VectorCross(VectorDerivative(angular_velocity(time), time), position_vector),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    angular_velocity_change_=angular_velocity,
    time_change_=time,
    position_vector=position_vector,
)
@validate_output(non_uniform_rotation_acceleration)
def calculate_non_uniform_rotation_acceleration(
    angular_velocity_change_: QuantityCoordinateVector,
    time_change_: Quantity,
    radius_vector_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    angular_velocity_ = (time / time_change_) * angular_velocity_change_

    result = law.rhs.subs(
        angular_velocity(time),
        angular_velocity_,
    ).doit().subs(
        position_vector,
        radius_vector_,
    )

    return QuantityCoordinateVector.from_expr(result)
