from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
)

# Description
# The Mach number is a dimensionless quantity in fluid dynamics representing
# the ratio of flow velocity past a boundary to the local speed of sound.
# Law: M = v / c, where
# v is velocity,
# c is speed of sound,
# M is Mach number.

velocity = Symbol("velocity", units.velocity)
speed_of_sound = Symbol("speed_of_sound", units.velocity)
mach_number = Symbol("mach_number", dimensionless)

law = Eq(mach_number, velocity / speed_of_sound)


def print_law() -> str:
    return print_expression(law)


@validate_input(velocity_=velocity, speed_of_sound_=speed_of_sound)
@validate_output(mach_number)
def calculate_mach_number(velocity_: Quantity, speed_of_sound_: Quantity) -> float:
    result_expr = solve(law, mach_number, dict=True)[0][mach_number]
    result_applied = result_expr.subs({
        velocity: velocity_,
        speed_of_sound: speed_of_sound_,
    })
    result = Quantity(result_applied)
    return convert_to_float(result)
