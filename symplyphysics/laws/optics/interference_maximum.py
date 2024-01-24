from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression,
                           validate_input, validate_output, dimensionless)

# Description
## If the stroke difference is equal to an integer number of wavelengths
## or an even number of half-waves, then a maximum intensity of interference will be observed.

# Law: delta = m * lambda
# Where:
## delta - optical travel difference for two coherent waves
## lambda - length of wave
## m - number of interference maximum


travel_difference = Symbol("travel_difference", units.length)
number_maximum = Symbol("number_maximum", dimensionless)
wave_length = Symbol("wave_length", units.length)

law = Eq(travel_difference, number_maximum * wave_length)


def print_law():
    return print_expression(law)


@validate_input(wave_length_=wave_length, number_maximum_=number_maximum)
@validate_output(travel_difference)
def calculate_travel_difference(wave_length_: Quantity, number_maximum_: Quantity | int) -> Quantity:
    solved = solve(law, wave_length, dict=True)[0][wave_length]
    result_expr = solved.subs({
        wave_length: wave_length_,
        number_maximum: number_maximum_
    })
    result = Quantity(result_expr)
    return result
