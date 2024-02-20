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
## The change in gas pressure when the temperature changes relative to 273 kelvin depends on temperature,
## standard pressure and thermal coefficient.
## https://www.fxyz.ru/формулы_по_физике/термодинамика_теория_теплоты/тепловое_расширение/расширение_газа/термический_коэффициент_давления/

## Law is: p = p0 * y * T, where
## p - pressure change,
## p0 - pressure at zero degrees celsius,
## y - thermal coefficient,
## T - temperature.

# Conditions:
## - the volume of gas remains constant when heated

pressure_change = Symbol("pressure", units.pressure)

standard_pressure = Symbol("standard_pressure", units.pressure)
thermal_coefficient = Symbol("thermal_coefficient", 1 / units.temperature)
temperature = Symbol("temperature", units.temperature)


law = Eq(pressure_change, standard_pressure * thermal_coefficient * temperature)


def print_law() -> str:
    return print_expression(law)


@validate_input(standard_pressure_=standard_pressure, thermal_coefficient_=thermal_coefficient, temperature_=temperature)
@validate_output(pressure_change)
def calculate_pressure_change(standard_pressure_: Quantity, thermal_coefficient_: Quantity, temperature_: Quantity) -> Quantity:
    result_expr = solve(law, pressure_change, dict=True)[0][pressure_change]
    result_expr = result_expr.subs({
        standard_pressure: standard_pressure_,
        thermal_coefficient: thermal_coefficient_,
        temperature: temperature_
    })
    return Quantity(result_expr)
