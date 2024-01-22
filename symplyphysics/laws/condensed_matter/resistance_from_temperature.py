from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The resistance depends on the temperature. For different materials, the value
## of the temperature coefficient and resistance at zero degrees celsius may differ.

## Law is: R = R0 * (1 + a * (T - 273.15)), where
## R - resistance,
## R0 - resistance at zero degrees celsius,
## a - temperature coefficient,
## T - temperature.

resistance = Symbol("resistance", units.impedance)

resistance_initial = Symbol("resistance_initial", units.impedance)
temperature_coefficient = Symbol("temperature_coefficient", 1 / units.temperature)
temperature = Symbol("temperature", units.temperature)

celsius_to_kelvin = Quantity(273.15 * units.kelvin)

law = Eq(resistance, resistance_initial * (1 + temperature_coefficient * (temperature - celsius_to_kelvin)))


def print_law() -> str:
    return print_expression(law)


@validate_input(resistance_initial_=resistance_initial,
    temperature_coefficient_=temperature_coefficient,
    temperature_=temperature)
@validate_output(resistance)
def calculate_resistance(resistance_initial_: Quantity, temperature_coefficient_: Quantity,
    temperature_: Quantity) -> Quantity:
    result_expr = solve(law, resistance, dict=True)[0][resistance]
    result_expr = result_expr.subs({
        resistance_initial: resistance_initial_,
        temperature_coefficient: temperature_coefficient_,
        temperature: temperature_
    })
    return Quantity(result_expr)
