from sympy import (Eq, solve)
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    quantities,
)

# Description
## The resistance depends on the temperature. For different materials, the value
## of the temperature coefficient and resistance at zero degrees celsius may differ.

## Law is: R = R0 * (1 + a * (T - T0)), where
## R - resistance,
## R0 - resistance at zero degrees celsius,
## a - temperature coefficient,
## T - temperature,
## T0 - 273.15 kelvin degrees.

# Links: BYJU's, similar formula for resistivity <https://byjus.com/physics/resistivity-temperature-dependence/>

resistance = Symbol("resistance", units.impedance)

resistance_initial = Symbol("resistance_initial", units.impedance)
temperature_coefficient = Symbol("temperature_coefficient", 1 / units.temperature)
temperature = symbols.temperature

law = Eq(
    resistance,
    resistance_initial * (1 + temperature_coefficient *
    (temperature - quantities.standard_conditions_temperature)))


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
