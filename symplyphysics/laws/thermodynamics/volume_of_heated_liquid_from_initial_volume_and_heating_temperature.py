from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
                           validate_output)

# Description
## As the temperature increases, the intensity of the thermal motion of the particles of matter increases. This leads to the fact that the molecules are more "actively" repelled from each other.
## As the empty spaces between molecules increase, most bodies increase in size when heated.

## Law: V = V_start * (1 + y * (t_end - t_start)
## Where:
## V_start is liquid volume at at the initial temperature
## V is volume of liquid at temperature t
## y is the coefficient of volumetric expansion of the liquid
## t_start is initial temperature in Kelvins
## t_end is final temperature in Kelvins


final_volume = Symbol("final_volume", units.volume)
start_volume = Symbol("start_volume", units.volume)
expansion_coefficient = Symbol("coefficient", 1 / units.temperature)
start_temperature = Symbol("start_temperature", units.temperature)
final_temperature = Symbol("final_temperature", units.temperature)


law = Eq(final_volume, start_volume * (1 + expansion_coefficient * (final_temperature - start_temperature)))


def print_law() -> str:
    return print_expression(law)


@validate_input(start_volume_=start_volume, expansion_coefficient_=expansion_coefficient, final_temperature_=final_temperature, start_temperature_=start_temperature)
@validate_output(final_volume)
def calculate_finish_volume(start_volume_, expansion_coefficient_, final_temperature_, start_temperature_: Quantity) -> Quantity:
    result_expr = solve(law, final_volume, dict=True)[0][final_volume]
    result_volume = result_expr.subs({
        start_volume: start_volume_,
        expansion_coefficient: expansion_coefficient_,
        final_temperature: final_temperature_,
        start_temperature: start_temperature_
    })
    return Quantity(result_volume)
