from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The Mach number is a quantity in fluid dynamics representing the ratio of flow velocity 
## past a boundary to the local speed of sound. The Mach number is primarily used to determine 
## the approximation with which a flow can be treated as an incompressible flow.

# Definition: M = u / c
## M - Mach number
## u - speed of moving object
## c - speed of sound in given medium

mach_number = Symbol("mach_number", dimensionless)
object_speed = Symbol("object_velocity", units.velocity)
speed_of_sound = Symbol("speed_of_sound", units.velocity)

definition = Eq(mach_number, object_speed / speed_of_sound)


def print_law() -> str:
    return print_expression(definition)


@validate_input(
    object_speed_=object_speed,
    speed_of_sound_=speed_of_sound,
)
@validate_output(mach_number)
def calculate_mach_number(
    object_speed_: Quantity,
    speed_of_sound_: Quantity,
) -> float:
    result = definition.rhs.subs({
        object_speed: object_speed_,
        speed_of_sound: speed_of_sound_,
    })
    return Quantity(result).scale_factor
