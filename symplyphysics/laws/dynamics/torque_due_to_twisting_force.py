from sympy import Eq, solve, sin
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    angle_type,
)

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


def print_law() -> str:
    return print_expression(law)


@validate_input(force_=force, distance_to_axis_=distance_to_axis, angle_=angle)
@validate_output(torque)
def calculate_torque(force_: Quantity, distance_to_axis_: Quantity, angle_: Quantity | float) -> Quantity:
    result = solve(law, torque)[0]
    angle_value = angle_.scale_factor if isinstance(angle_, Quantity) else angle_
    result_torque = result.subs({
        distance_to_axis: distance_to_axis_,
        force: force_,
        angle: angle_value,
    })
    return Quantity(result_torque)
