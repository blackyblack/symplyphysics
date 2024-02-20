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

body_velocity_in_first_frame = Symbol("body_velocity_in_first_frame", units.velocity)
body_velocity_in_second_frame = Symbol("body_velocity_in_second_frame", units.velocity)
second_frame_velocity_in_first_frame = Symbol(
    "second_frame_velocity_in_first_frame", units.velocity
)

law = Eq(
    body_velocity_in_first_frame,
    body_velocity_in_second_frame + second_frame_velocity_in_first_frame,
)

def print_law() -> str:
    return print_expression(law)


@validate_input(
    velocity_in_second_frame_=body_velocity_in_second_frame,
    second_frame_velocity_in_first_frame_=second_frame_velocity_in_first_frame,
)
@validate_output(body_velocity_in_first_frame)
def calculate_body_velocity_in_first_frame(
    velocity_in_second_frame_: Quantity,
    second_frame_velocity_in_first_frame_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        body_velocity_in_second_frame: velocity_in_second_frame_,
        second_frame_velocity_in_first_frame: second_frame_velocity_in_first_frame_,
    })
    return Quantity(result)
