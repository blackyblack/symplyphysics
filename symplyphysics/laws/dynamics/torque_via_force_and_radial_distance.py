"""
Torque via force and radial distance
====================================

*Torque* is a turning action on a body about a rotation axis due to a force.

**Notes:**

#. The position vector of a point in space, also known as location or radius vector,
   is the vector connecting the origin of the coordinate system and the given point.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Torque#Definition_and_relation_to_other_physical_quantities>`__.
"""

from sympy import Eq, solve, sin, symbols as sympy_symbols
from symplyphysics import symbols, Quantity, validate_input, validate_output
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.symbols.quantities import scale_factor
from symplyphysics.laws.dynamics.vector import (
    torque_vector_of_twisting_force as _torque_vector_def,)
from symplyphysics.laws.geometry.vector import (
    dot_product_is_proportional_to_cosine_between_vectors as _dot_product_law,)

from symplyphysics.core.experimental.solvers import solve_for_vector
from symplyphysics.core.experimental.coordinate_systems import CoordinateVector, CARTESIAN
from symplyphysics.core.experimental.vectors import VectorNorm

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

_force_vector = CoordinateVector(sympy_symbols("F_x:z", real=True), CARTESIAN)
_position_vector = CoordinateVector(sympy_symbols("x:z", real=True), CARTESIAN)

_torque_vector_derived = solve_for_vector(
    _torque_vector_def.law,
    _torque_vector_def.torque,
).subs({
    _torque_vector_def.force: _force_vector,
    _torque_vector_def.position_vector: _position_vector,
})

_torque_magnitude_derived = VectorNorm(_torque_vector_derived)

_angle_between_vectors = solve(
    _dot_product_law.law,
    _dot_product_law.angle_between_vectors,
)[-1].subs({
    _dot_product_law.first_vector: _force_vector,
    _dot_product_law.second_vector: _position_vector,
})

_torque_magnitude_expected = solve(law, torque)[0].subs({
    force: VectorNorm(_force_vector),
    radial_distance: VectorNorm(_position_vector),
    angle_between_vectors: _angle_between_vectors,
})

assert expr_equals(_torque_magnitude_derived, _torque_magnitude_expected)


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


# UNIQUE_LAW_ID: 214
