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

# Description
## A body is rotating about a fixed axis with constant angular acceleration. Its angular
## displacement from initial position is a quadratic function of time and depends on
## its initial angular velocity and angular acceleration.

# Law: theta = w0 * t + alpha * t**2 / 2
## theta - angular displacement from initial position
## w0 - initial angular velocity
## alpha - constant angular acceleration
## t - time

## Conditions:
## Angular acceleration of the body is constant.

angular_displacement = Symbol("angular_displacement", angle_type)
initial_angular_velocity = Symbol("initial_angular_velocity", angle_type / units.time)
angular_acceleration = Symbol("angular_acceleration", angle_type / units.time**2)
time = Symbol("time", units.time)

law = Eq(
    angular_displacement,
    initial_angular_velocity * time + angular_acceleration * time**2 / 2,
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    initial_angular_velocity_=initial_angular_velocity,
    angular_acceleration_=angular_acceleration,
    time_=time,
)
@validate_output(angular_displacement)
def calculate_angular_displacement(
    initial_angular_velocity_: Quantity,
    angular_acceleration_: Quantity,
    time_: Quantity,
) -> Quantity | float:
    result = law.rhs.subs({
        initial_angular_velocity: initial_angular_velocity_,
        angular_acceleration: angular_acceleration_,
        time: time_,
    })
    return Quantity(result)
