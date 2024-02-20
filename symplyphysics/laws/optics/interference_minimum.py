from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless)

# Description
## If the travel difference is equal to an odd number of half-waves,
## then at this point of the screen there will be a minimum of intensity during interference.

# Law: delta = (2m+1) * lambda / 2
# Where:
## delta - optical travel difference for two coherent waves
## lambda - length of waves (wave 1 and wave 2)
## m - number of interference minimum (m = 0, 1, 2, ...)

# Condition: The two waves must be coherent

travel_difference = Symbol("travel_difference", units.length)
number_minimum = Symbol("number_minimum", dimensionless)
wave_length = Symbol("wave_length", units.length)

law = Eq(travel_difference, (2 * number_minimum + 1) * wave_length / 2)


def print_law() -> str:
    return print_expression(law)


@validate_input(wave_length_=wave_length, number_minimum_=number_minimum)
@validate_output(travel_difference)
def calculate_travel_difference(wave_length_: Quantity, number_minimum_: int) -> Quantity:
    solved = solve(law, travel_difference, dict=True)[0][travel_difference]
    result_expr = solved.subs({wave_length: wave_length_, number_minimum: number_minimum_})
    result = Quantity(result_expr)
    return result
