from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, dimensionless)

# Description
## If the difference in the course of two waves is equal to an integer number of waves
## (i.e. an even number of half-waves) Δ = mλ, where m = 0, 1, 2, ...,
## then an interference maximum is formed at the point of superposition of these waves.

# Links:
## https://vogueindustry.com/17289645-interference-patterns-maximum-and-minimum-conditions#menu-9
## https://www.livelaptopspec.com/what-is-maximum-constructive-interference/

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


@validate_input(wave_length_=wave_length, number_maximum_=number_maximum)
@validate_output(travel_difference)
def calculate_travel_difference(wave_length_: Quantity, number_maximum_: int) -> Quantity:
    solved = solve(law, travel_difference, dict=True)[0][travel_difference]
    result_expr = solved.subs({wave_length: wave_length_, number_maximum: number_maximum_})
    return Quantity(result_expr)
