from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, dimensionless, validate_input,
    validate_output)

# Description
## There is an empirical formula for the spectral lines of the hydrogen atom.

# Law: f = R * (1 / m^2 - 1 / n^2 ), where
## Vmedium is speed of electromagnetic wave in medium,
## f - frequency of the electron transition in hydrogen,
## R - Rydberg's constant,
## m - the number of the level to which the electron goes,
## n - the number of the level from which the electron leaves.

transition_frequency = Symbol("transition_frequency", units.frequency)

number_level_to  = Symbol("number_level_to", dimensionless)
number_level_from = Symbol("number_level_from", dimensionless)

rydberg_constant = Quantity(3.29e15 * units.hertz)

law = Eq(transition_frequency, rydberg_constant * ((1 / number_level_to**2) - (1 / number_level_from**2)))


def print_law() -> str:
    return print_expression(law)


@validate_input(number_level_to_=number_level_to, number_level_from_=number_level_from)
@validate_output(transition_frequency)
def calculate_transition_frequency(number_level_to_: float, number_level_from_: float) -> Quantity:
    result_expr = solve(law, transition_frequency, dict=True)[0][transition_frequency]
    transition_frequency_applied = result_expr.subs({
        number_level_to: number_level_to_,
        number_level_from: number_level_from_
    })
    return Quantity(transition_frequency_applied)
