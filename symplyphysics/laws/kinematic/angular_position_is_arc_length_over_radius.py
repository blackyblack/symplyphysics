from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    angle_type,
)

# Description:
## To describe the rotation of a rigid body about a fixed axis (rotational axis), a reference line is
## assumed to be fixed in the body, perpendicular to that axis and rotating with the body. The angular
## position of this is measured relative to a fixed direction and is expressed as a ratio of the arc 
## length of a circular path and its radius r.

# Law: theta = s / r
## theta - angular position
## s - arc length
## r - path radius

angular_position = Symbol("angular_position", angle_type)
arc_length = Symbol("arc_length", units.length)
path_radius = Symbol("path_radius", units.length)

law = Eq(angular_position, arc_length / path_radius)


def print_law() -> str:
    return print_expression(law)


@validate_input(arc_length_=arc_length, path_radius_=path_radius)
@validate_output(angular_position)
def calculate_angular_position(arc_length_: Quantity, path_radius_: Quantity) -> Quantity:
    result = law.rhs.subs({
        arc_length: arc_length_,
        path_radius: path_radius_,
    })
    result_angle = Quantity(result)
    return result_angle
