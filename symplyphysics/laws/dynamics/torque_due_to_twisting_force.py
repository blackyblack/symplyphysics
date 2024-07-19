"""
Torque due to twisting force
============================

*Torque* is a turning action on a body about a rotation axis due to a force.

**Notes:**

#. The position vector of a point in space, also known as location or radius vector,
   is the vector connecting the origin of the coordinate system and the given point.
"""

from sympy import Eq, solve, sin, symbols as sympy_symbols, sqrt
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    angle_type,
    Vector,
    dot_vectors,
    vector_magnitude,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.symbols.quantities import scale_factor
from symplyphysics.laws.dynamics.vector import torque_vector_of_twisting_force as torque_vector_def

torque = Symbol("torque", units.force * units.length)
r"""
The magnitude of the torque applied at the given point.

Symbol:
    tau

Latex:
    :math:`\tau`
"""

force = symbols.dynamics.force
"""
The magnitude of the :attr:`~symplyphysics.symbols.dynamics.force` exerted on the given point.

Symbol:
    F
"""

distance_to_axis = Symbol("distance_to_axis", units.length)
"""
The distance to axis from the given point.

Symbol:
    r
"""

angle = Symbol("angle", angle_type)
r"""
The angle between the position vector of the given point and the force vector.

Symbol:
    phi

Latex:
    :math:`\varphi`
"""

law = Eq(torque, distance_to_axis * force * sin(angle))
r"""
tau = F * r * sin(phi)

Latex:
    :math:`\tau = F r \sin{\varphi}`
"""

# Derive law from its vector counterpart.
force_vector = Vector(sympy_symbols("force_x:z", real=True))
position_vector = Vector(sympy_symbols("x:z", real=True))

torque_vector_derived = torque_vector_def.torque_definition(force_vector, position_vector)
torque_magnitude_derived = vector_magnitude(torque_vector_derived)

force_magnitude = vector_magnitude(force_vector)
position_magnitude = vector_magnitude(position_vector)

# Use the definition of dot product (a, b) = |a| * |b| * cos(a, b) to find the sine of angle between vectors
cosine_of_angle_in_between = (dot_vectors(force_vector, position_vector) / force_magnitude /
    position_magnitude)
sine_of_angle_in_between = sqrt(1 - cosine_of_angle_in_between**2)

torque_magnitude_from_law = solve(law, torque)[0].subs({
    force: force_magnitude,
    distance_to_axis: position_magnitude,
    sin(angle): sine_of_angle_in_between,
})

assert expr_equals(torque_magnitude_derived, torque_magnitude_from_law)


@validate_input(force_=force, distance_to_axis_=distance_to_axis, angle_=angle)
@validate_output(torque)
def calculate_torque(force_: Quantity, distance_to_axis_: Quantity,
    angle_: Quantity | float) -> Quantity:
    result = solve(law, torque)[0]
    angle_value = scale_factor(angle_)
    result_torque = result.subs({
        distance_to_axis: distance_to_axis_,
        force: force_,
        angle: angle_value,
    })
    return Quantity(result_torque)
