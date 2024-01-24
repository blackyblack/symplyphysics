from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression,
                           validate_input, validate_output, dimensionless)

# Description
## If the stroke difference is equal to an odd number of half-waves,
## then at this point of the screen there will be a minimum of intensity during interference.

# Law: delta = (2m+1) * lambda / 2
# Where:
## delta - optical travel difference for two coherent waves
## lambda - length of wave
## m - number of interference minimum


travel_difference = Symbol("travel_difference", units.length)
number_minimum = Symbol("number_minimum", dimensionless)
wave_length = Symbol("wave_length", units.length)

law = Eq(travel_difference, number_minimum * wave_length)


def print_law():
    return print_expression(law)


@validate_input(wave_length_=wave_length, number_minimum_=number_minimum)
@validate_output(travel_difference)
def calculate_travel_difference(wave_length_: Quantity, number_minimum_: Quantity | int) -> Quantity:
    solved = solve(law, wave_length, dict=True)[0][wave_length]
    result_expr = solved.subs({
        wave_length: wave_length_,
        number_minimum: number_minimum_
    })
    result = Quantity(result_expr)
    return result