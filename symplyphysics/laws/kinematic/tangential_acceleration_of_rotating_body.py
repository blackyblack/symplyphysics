from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    angle_type
)

# Description
## The tangential acceleration of a rotating body represents the change in magnitude
## of the velocity vector, and is tangent to the path of the body. It is proportional
## to the angular acceleration of the body and its rotation radius.

# Law: a_t = alpha * r
## a_t - tangential acceleration
## alpha - angular acceleration
## r - rotation radius

tangential_acceleration = Symbol("tangential_acceleration", units.acceleration)
angular_acceleration = Symbol("angular_acceleration", angle_type / units.time**2)
rotation_radius = Symbol("rotation_radius", units.length)

law = Eq(tangential_acceleration, angular_acceleration * rotation_radius)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    angular_acceleration_=angular_acceleration,
    rotation_radius_=rotation_radius,
)
@validate_output(tangential_acceleration)
def calculate_tangential_acceleration(
    angular_acceleration_: Quantity, 
    rotation_radius_: Quantity
) -> Quantity:
    result_expr = solve(law, tangential_acceleration)[0]
    result = result_expr.subs({
        angular_acceleration: angular_acceleration_,
        rotation_radius: rotation_radius_,
    })
    return Quantity(result)
