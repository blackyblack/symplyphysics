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
## If a body is rotating about a fixed axis with constant angular velocity, its angular position
## is linearly proportional to time and angular velocity.

# Law: theta = theta0 + w * t
## theta - end angle
## theta0 - start angle
## w - constant angular velocity
## t - time

final_angular_position = Symbol("final_angular_position", angle_type)
initial_angular_position = Symbol("initial_angular_position", angle_type)
angular_velocity = Symbol("angular_velocity", angle_type / units.time)
time = Symbol("time", units.time)

law = Eq(
    final_angular_position,
    initial_angular_position + angular_velocity * time,
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    initial_angular_position_=initial_angular_position,
    angular_velocity_=angular_velocity,
    time_=time,
)
@validate_output(final_angular_position)
def calculate_angular_position(
    initial_angular_position_: Quantity | float,
    angular_velocity_: Quantity,
    time_: Quantity,
) -> Quantity:
    initial_angular_position_ = (
        initial_angular_position_.scale_factor
        if isinstance(initial_angular_position_, Quantity)
        else initial_angular_position_
    )
    result = law.rhs.subs({
        initial_angular_position: initial_angular_position_,
        angular_velocity: angular_velocity_,
        time: time_,
    })
    return Quantity(result)
