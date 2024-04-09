from sympy import Eq, solve, pi, Abs
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    angle_type,
)
from symplyphysics.core.convert import convert_to_dimensionless
# Description
## The phase difference of the oscillations of two points spaced at distances r1 and r2 from the source of
## the oscillations depends on these distances and wavelength.

## Law is: phi = 2 * pi * |r2 - r1| / L, where
## phi - phase difference,
## r2 - distance to the second point,
## r1 - distance to the first point,
## L - wavelength,
## || - absolute value.

phase_difference = Symbol("phase_difference", angle_type)

distance_first_point = Symbol("distance_first_point", units.length)
distance_second_point = Symbol("distance_second_point", units.length)
wavelength = Symbol("wavelength", units.length)

law = Eq(phase_difference, 2 * pi * Abs(distance_second_point - distance_first_point) / wavelength)


def print_law() -> str:
    return print_expression(law)


@validate_input(distance_first_point_=distance_first_point,
    distance_second_point_=distance_second_point,
    wavelength_=wavelength)
@validate_output(phase_difference)
def calculate_phase_difference(distance_first_point_: Quantity, distance_second_point_: Quantity,
    wavelength_: Quantity) -> float:
    result_expr = solve(law, phase_difference, dict=True)[0][phase_difference]
    result_expr = result_expr.subs({
        distance_first_point: distance_first_point_,
        distance_second_point: distance_second_point_,
        wavelength: wavelength_
    }).doit()
    return convert_to_dimensionless(Quantity(result_expr))
