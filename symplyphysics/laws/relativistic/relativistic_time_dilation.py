from sympy import Eq, solve, sqrt
from sympy.physics.units import speed_of_light

from symplyphysics import (Quantity, Symbol, print_expression, units, validate_input,
                           validate_output)


# Description
# Time dilation is the difference in elapsed time as measured by two clocks,
# either due to a relative velocity between them
# Law: t_rel = t * sqrt(1 - v**2 / c**2), where
# t_rel is time elapsed for moving observer,
# t is is time elapsed for rest observer,
# v is velocity,
# c is speed of light.

rest_time = Symbol("rest_time", units.time)
velocity = Symbol("velocity", units.velocity)
relativistic_time = Symbol("relativistic_time", units.time)

law = Eq(relativistic_time, rest_time /
         sqrt(1 - velocity**2 / speed_of_light**2))


def print_law() -> str:
    return print_expression(law)


@validate_input(rest_time_=rest_time, velocity_=velocity)
@validate_output(relativistic_time)
def calculate_relativistic_time(rest_time_: Quantity, velocity_: Quantity) -> Quantity:
    result_expr = solve(law, relativistic_time, dict=True)[
        0][relativistic_time]
    time_applied = result_expr.subs(
        {rest_time: rest_time_, velocity: velocity_})
    return Quantity(time_applied)
