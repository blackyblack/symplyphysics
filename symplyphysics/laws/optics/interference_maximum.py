from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression,
                           validate_input, validate_output, dimensionless)

# Description
## If the stroke difference is equal to an even number of half-waves,
## then the interference intensity will be maximum at this point of the screen.

# Law: delta = m * lambda
# Where:
## delta - optical travel difference for two coherent waves
## lambda - length of waves (wave 1 and wave 2)
## m - number of interference maximum (m = 0, 1, 2, ...)

# Condition: The two waves must be coherent


travel_difference = Symbol("travel_difference", units.length)
number_maximum = Symbol("number_maximum", dimensionless)
wave_length = Symbol("wave_length", units.length)

law = Eq(travel_difference, number_maximum * wave_length)


def print_law():
    return print_expression(law)


@validate_input(wave_length_=wave_length, number_maximum_=number_maximum)
@validate_output(travel_difference)
def calculate_travel_difference(wave_length_: Quantity, number_maximum_: int) -> Quantity:
    solved = solve(law, travel_difference, dict=True)[0][travel_difference]
    result_expr = solved.subs({
        wave_length: wave_length_,
        number_maximum: number_maximum_
    })
    result = Quantity(result_expr)
    return result
