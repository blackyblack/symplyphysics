"""
Centripetal acceleration via vector rejection
=============================================

*Centripetal acceleration* is the acceleration of a body in a rotating coordinate system
which is directed towards the axis of rotation.

Also see :doc:`laws.kinematics.vector.centripetal_acceleration_via_cross_product`.

**Notation:**

#. :math:`\\left( \\vec a, \\vec b \\right)` (:code:`dot(a, b)`) is the dot product between
   vectors :math:`\\vec a` and :math:`\\vec b`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Centripetal_force#Derivation_using_vectors>`__.
"""

from sympy import Eq
from symplyphysics import validate_input, validate_output, symbols

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorDot
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

centripetal_acceleration = clone_as_vector_symbol(
    symbols.acceleration,
    display_symbol="a_cp",
    display_latex="{\\vec a}_\\text{cp}",
)
"""
Vector of the body's centripetal :symbols:`acceleration`.
"""

angular_velocity = clone_as_vector_symbol(symbols.angular_speed)
"""
Pseudovector of the angular velocity of the body's rotation. See :symbols:`angular_speed`.
"""

position_vector = clone_as_vector_symbol(symbols.distance_to_origin)
"""
The body's position vector. See :symbols:`distance_to_origin`.
"""

# NOTE: the following is an equivalent form of the original `a_cp = -1 * dot(w, w) * reject(r, w)`
# where `reject(r, w) = r - project(r, w)` and `project(r, w) = w * dot(r, w) / dot(w, w)`.
# Maybe implement `VectorReject` and `VectorProject` classes if enough laws need them.
law = Eq(
    centripetal_acceleration,
    angular_velocity * VectorDot(position_vector, angular_velocity) -
    position_vector * VectorDot(angular_velocity, angular_velocity),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    angular_velocity_=angular_velocity,
    radius_vector_=position_vector,
)
@validate_output(centripetal_acceleration)
def calculate_centripetal_acceleration(
    angular_velocity_: QuantityCoordinateVector,
    radius_vector_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = law.rhs.subs({
        angular_velocity: angular_velocity_,
        position_vector: radius_vector_,
    })

    return QuantityCoordinateVector.from_expr(result)
