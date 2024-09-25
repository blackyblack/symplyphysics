"""
Torque via force and radial distance
====================================

*Torque* is a turning action on a body about a rotation axis due to a force.

**Notes:**

#. The position vector of a point in space, also known as location or radius vector,
   is the vector connecting the origin of the coordinate system and the given point.
"""

from sympy import Eq, solve, sin, symbols as sympy_symbols, sqrt
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
    Vector,
    dot_vectors,
    vector_magnitude,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.symbols.quantities import scale_factor
from symplyphysics.laws.dynamics.vector import torque_vector_of_twisting_force as torque_vector_def

torque = symbols.torque
"""
The magnitude of the :symbols:`torque` applied at the given point.
"""

force = symbols.force
"""
The magnitude of the :symbols:`force` exerted on the given point.
"""

radial_distance = symbols.distance_to_axis
"""
The :symbols:`distance_to_axis` from the given point.
"""

angle_between_vectors = symbols.angle
"""
The :symbols:`angle` between the position vector of the given point and the force vector.
"""

law = Eq(torque, radial_distance * force * sin(angle_between_vectors))
"""
:laws:symbol::

:laws:latex::
"""

# Derive law from its vector counterpart.
_force_vector = Vector(sympy_symbols("force_x:z", real=True))
_position_vector = Vector(sympy_symbols("x:z", real=True))

_torque_vector_derived = torque_vector_def.torque_definition(_force_vector, _position_vector)
_torque_magnitude_derived = vector_magnitude(_torque_vector_derived)

_force_magnitude = vector_magnitude(_force_vector)
_position_magnitude = vector_magnitude(_position_vector)

# Use the definition of dot product (a, b) = |a| * |b| * cos(a, b) to find the sine of angle_between_vectors between vectors
_cosine_of_angle_in_between = (dot_vectors(_force_vector, _position_vector) / _force_magnitude /
    _position_magnitude)
_sine_of_angle_in_between = sqrt(1 - _cosine_of_angle_in_between**2)

_torque_magnitude_from_law = solve(law, torque)[0].subs({
    force: _force_magnitude,
    radial_distance: _position_magnitude,
    sin(angle_between_vectors): _sine_of_angle_in_between,
})

assert expr_equals(_torque_magnitude_derived, _torque_magnitude_from_law)


@validate_input(force_=force, distance_to_axis_=radial_distance, angle_=angle_between_vectors)
@validate_output(torque)
def calculate_torque(force_: Quantity, distance_to_axis_: Quantity,
    angle_: Quantity | float) -> Quantity:
    result = solve(law, torque)[0]
    angle_value = scale_factor(angle_)
    result_torque = result.subs({
        radial_distance: distance_to_axis_,
        force: force_,
        angle_between_vectors: angle_value,
    })
    return Quantity(result_torque)
