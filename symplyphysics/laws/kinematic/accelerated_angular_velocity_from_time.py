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
## If a body is rotating about a fixed axis with constant angular acceleration, its angular
## velocity is a linear function of time.

# Law: w = w0 + alpha * t
## w - angular velocity at time t
## w0 - initial angular velocity
## alpha - constant angular acceleration
## t - time

# Conditions
## Angular acceleration should be constant.

angular_velocity = Symbol("angular_velocity", angle_type / units.time)
initial_angular_velocity = Symbol("initial_angular_velocity", angle_type / units.time)
angular_acceleration = Symbol("angular_acceleration", angle_type / units.time**2)
time = Symbol("time", units.time)

law = Eq(angular_velocity, initial_angular_velocity + angular_acceleration * time)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    initial_angular_velocity_=initial_angular_velocity,
    angular_acceleration_=angular_acceleration,
    time_=time,
)
@validate_output(angular_velocity)
def calculate_angular_velocity(
    initial_angular_velocity_: Quantity,
    angular_acceleration_: Quantity,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        initial_angular_velocity: initial_angular_velocity_,
        angular_acceleration: angular_acceleration_,
        time: time_,
    })
    return Quantity(result)
