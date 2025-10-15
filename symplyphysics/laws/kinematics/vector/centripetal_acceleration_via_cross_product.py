"""
Centripetal acceleration via cross product
==========================================

*Centripetal acceleration* is the acceleration of a body in a rotating coordinate system which is
directed towards the axis of rotation.

Also see :doc:`laws.kinematics.vector.centrifugal_acceleration_via_centripetal_acceleration`.

**Notation:**

#. :math:`\\left[ \\vec a, \\vec b \\right]` (:code:`cross(a, b)`) is the cross product between
   :math:`\\vec a` and :math:`\\vec b`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Centripetal_force#Derivation_using_vectors>`__.
"""

from sympy import Eq, symbols as sym_symbols
from symplyphysics import validate_input, validate_output, symbols

from symplyphysics.laws.kinematics.vector import (
    centripetal_acceleration_via_vector_rejection as rejection_law,)

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorCross
from symplyphysics.core.experimental.coordinate_systems import (CoordinateVector, CARTESIAN,
    QuantityCoordinateVector)
from symplyphysics.core.experimental.solvers import vector_equals

centripetal_acceleration = clone_as_vector_symbol(
    symbols.acceleration,
    display_symbol="a_cp",
    display_latex="{\\vec a}_\\text{cp}",
)
"""
Vector of centripetal :symbols:`acceleration`.
"""

angular_velocity = clone_as_vector_symbol(symbols.angular_speed)
"""
Pseudovector of angular velocity of the body. See :symbols:`angular_speed`.
"""

position_vector = clone_as_vector_symbol(symbols.distance_to_origin)
"""
Position vector of the body. See :symbols:`distance_to_origin`.
"""

law = Eq(
    centripetal_acceleration,
    VectorCross(angular_velocity, VectorCross(angular_velocity, position_vector)),
)
"""
:laws:symbol::

:laws:latex::
"""

# Prove that obtaining centripetal acceleration via cross product and via vector rejection yield
# same result.

_angular_velocity = CoordinateVector(sym_symbols("angular_velocity_x:z"), CARTESIAN)
_position_vector = CoordinateVector(sym_symbols("position_vector_x:z"), CARTESIAN)

_cross_product_result = law.rhs.subs({
    angular_velocity: _angular_velocity,
    position_vector: _position_vector,
})
_cross_product_result = CoordinateVector.from_expr(_cross_product_result)

_rejection_result = rejection_law.law.rhs.subs({
    rejection_law.angular_velocity: _angular_velocity,
    rejection_law.position_vector: _position_vector,
})

_rejection_result = CoordinateVector.from_expr(_rejection_result)

assert vector_equals(_cross_product_result, _rejection_result)


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


# UNIQUE_LAW_ID: 447
