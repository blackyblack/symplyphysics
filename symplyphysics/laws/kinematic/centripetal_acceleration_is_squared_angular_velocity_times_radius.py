from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    angle_type,
    validate_input,
    validate_output,
)

# Description
## Centripetal acceleration is defined as the change in velocity tangential to the velocity vector.
## See more at [its definition through linear velocity](./centripetal_acceleration_is_squared_velocity_by_radius.py)

# Law: a_n = w**2 * r
## a_n - centripetal (or normal) acceleration
## w - angular velocity
## r - curve radius

centripetal_acceleration = Symbol("centripetal_acceleration", units.acceleration)
angular_velocity = Symbol("angular_velocity", angle_type / units.time)
curve_radius = Symbol("curve_radius", units.length)

law = Eq(centripetal_acceleration, angular_velocity**2 * curve_radius)


def print_law() -> str:
    return print_expression(law)


@validate_input(angular_velocity_=angular_velocity, curve_radius_=curve_radius)
@validate_output(centripetal_acceleration)
def calculate_centripetal_acceleration(angular_velocity_: Quantity, curve_radius_: Quantity) -> Quantity:
    result = law.rhs.subs({
        angular_velocity: angular_velocity_,
        curve_radius: curve_radius_,
    })
    return Quantity(result)
