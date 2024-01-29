from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    angle_type,
    Vector,
    vector_magnitude,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematic.vector import (
    linear_displacement_is_angular_displacement_cross_radius as displacement_law
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


# Derive law from its [vector counterpart](../linear_displacement_is_angular_displacement_cross_radius.py)

angle_x = Symbol("angle_x", angle_type)
angle_y = Symbol("angle_y", angle_type)
angle_z = Symbol("angle_z", angle_type)
angle_vec = Vector([angle_x, angle_y, angle_z])
angle_norm = vector_magnitude(angle_vec)

radius_x = Symbol("radius_x", units.length)
radius_y = Symbol("radius_y", units.length)
radius_z = Symbol("radius_z", units.length)
radius_vec = Vector([radius_x, radius_y, radius_z])
radius_norm = vector_magnitude(radius_vec)

displacement_from_law = solve(
    law, arc_length
)[0].subs({
    angular_position: angle_norm,
    path_radius: radius_norm,
})

# From [the law](../linear_displacement_is_angular_displacement_cross_radius.py), taking the vector 
# norm of both sides, we have norm(s) = norm(theta) * norm(r) * sin(phi) where s is the linear 
# displacement vector, theta is the angular displacement pseudovector, r is the radius vector, and 
# phi is the and between the last two. As a necessary condition, theta and r are orthogonal to each 
# other, therefore phi is pi/2 and sin(phi) is 1, and so abs(s) = abs(theta) * abs(r).

displacement_derived = angle_norm * radius_norm

assert expr_equals(displacement_from_law, displacement_derived)


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
