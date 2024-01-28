from sympy import Eq, solve, sin, symbols, sqrt
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    angle_type,
    Vector,
    dot_vectors,
    vector_magnitude,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics.vector import torque_vector_of_twisting_force as torque_vector_def

# Description
## Torque is a turning action on a body about a rotation axis due to a force.

# Law: tau = F * r * sin(phi)
## tau - torque
## F - force exerted at the given point
## r - distance to axis from the given point
## phi - angle between the position vector of the given point and force vector

# Note
## The position vector of a point in space, also known as location or radius vector,
## is the vector connecting the origin of the coordinate system and the given point.

torque = Symbol("torque", units.force * units.length)
force = Symbol("force", units.force)
distance_to_axis = Symbol("distance_to_axis", units.length)
angle = Symbol("angle", angle_type)

law = Eq(torque, distance_to_axis * force * sin(angle))

# Derive law from its vector counterpart.
force_vector = Vector(symbols("force_x:z", real=True))
position_vector = Vector(symbols("x:z", real=True))

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


def print_law() -> str:
    return print_expression(law)


@validate_input(force_=force, distance_to_axis_=distance_to_axis, angle_=angle)
@validate_output(torque)
def calculate_torque(force_: Quantity, distance_to_axis_: Quantity,
    angle_: Quantity | float) -> Quantity:
    result = solve(law, torque)[0]
    angle_value = angle_.scale_factor if isinstance(angle_, Quantity) else angle_
    result_torque = result.subs({
        distance_to_axis: distance_to_axis_,
        force: force_,
        angle: angle_value,
    })
    return Quantity(result_torque)
