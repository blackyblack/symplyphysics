from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The law of classical addition of velocities, usually attributed to Galileo and called the
## Galilean law of velocity addition, states that the velocity of a body in an inertial reference
## frame A can be found as a sum of its velocity in another inertial reference frame B and the
## velocity of frame B relative to frame A.

# Law: u = v + u'
## u - velocity of body relative to frame A
## u' - velocity of body relative to frame B
## v - velocity of frame B relative to frame A

# Conditions
## - Space and time are absolute
## - Applies to inertial frames of reference

speed_relative_to_first_frame = Symbol("velocity_relative_to_first_frame", units.velocity)
speed_relative_to_second_frame = Symbol("velocity_relative_to_second_frame", units.velocity)
relative_speed_of_second_frame_relative_to_first = Symbol(
    "relative_speed_of_second_frame_relative_to_first", units.velocity
)

law = Eq(
    speed_relative_to_first_frame,
    speed_relative_to_second_frame + relative_speed_of_second_frame_relative_to_first,
)

def print_law() -> str:
    return print_expression(law)


@validate_input(
    speed_relative_to_second_frame_=speed_relative_to_second_frame,
    relative_speed_of_second_frame_relative_to_first_=relative_speed_of_second_frame_relative_to_first,
)
@validate_output(speed_relative_to_first_frame)
def calculate_speed_relative_to_first_frame(
    speed_relative_to_second_frame_: Quantity,
    relative_speed_of_second_frame_relative_to_first_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        speed_relative_to_second_frame: speed_relative_to_second_frame_,
        relative_speed_of_second_frame_relative_to_first: relative_speed_of_second_frame_relative_to_first_,
    })
    return Quantity(result)
