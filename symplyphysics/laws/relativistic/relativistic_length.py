from sympy import Eq, solve, sqrt
from sympy.physics.units import speed_of_light

from symplyphysics import (Quantity, Symbol, print_expression, units,
                           validate_input, validate_output)

# Description
# Length contraction is the phenomenon that a moving object's length is
# measured to be shorter than its proper length, which is the length as
# measured in the object's own rest frame.
# Law: l_rel = l / sqrt(1 - v**2 / c**2), where
# l is rest length,
# v is velocity,
# c is speed of light,
# l_rel is relativistic length.


rest_length = Symbol("rest_length", units.length)
velocity = Symbol("velocity", units.velocity)
relativistic_length = Symbol("relativistic_length", units.length)

law = Eq(relativistic_length, rest_length *
         sqrt(1 - velocity**2 / speed_of_light**2))


def print_law():
    print_expression(law)


@validate_input(rest_length_=rest_length, velocity_=velocity)
@validate_output(relativistic_length)
def calculate_relativistic_length(rest_length_: Quantity, velocity_: Quantity) -> Quantity:
    result_expr = solve(law, relativistic_length, dict=True)[
        0][relativistic_length]
    mass_applied = result_expr.subs(
        {rest_length: rest_length_, velocity: velocity_})
    return Quantity(mass_applied)
