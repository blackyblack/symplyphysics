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
## The gas standard_pressure at a certain temperature depends on the temperature, standard standard_pressure and thermal coefficient.

## Law is: p = p0 * y * T, where
## p - pressure,
## p0 - pressure at zero degrees celsius,
## a - thermal coefficient,
## T - temperature.

pressure = Symbol("pressure", units.pressure)

standard_pressure = Symbol("standard_pressure", units.pressure)
thermal_coefficient = Symbol("thermal_coefficient", 1 / units.temperature)
temperature = Symbol("temperature", units.temperature)


law = Eq(pressure, standard_pressure * thermal_coefficient * temperature)


def print_law() -> str:
    return print_expression(law)


@validate_input(standard_pressure_=standard_pressure, thermal_coefficient_=thermal_coefficient, temperature_=temperature)
@validate_output(pressure)
def calculate_pressure(standard_pressure_: Quantity, thermal_coefficient_: Quantity, temperature_: Quantity) -> Quantity:
    result_expr = solve(law, pressure, dict=True)[0][pressure]
    result_expr = result_expr.subs({
        standard_pressure: standard_pressure_,
        thermal_coefficient: thermal_coefficient_,
        temperature: temperature_
    })
    return Quantity(result_expr)
