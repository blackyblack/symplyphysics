from sympy import (Eq, solve, S)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
    convert_to
)

# Description
## Relative humidity is the ratio of water vapor pressure p to saturated vapor pressure p0 at a given temperature.

## Law is: ph = (p / p_sat) * 100, where
## ph - relative humidity,
## p - water vapor pressure,
## p_sat - saturated vapor pressure.

# Conditions:
## - the process is isobaric-isothermal.

relative_humidity = Symbol("relative_humidity", dimensionless)

water_vapor_pressure = Symbol("water_vapor_pressure", units.pressure)
saturated_vapor_pressure = Symbol("saturated_vapor_pressure", units.pressure)

law = Eq(relative_humidity, (water_vapor_pressure / saturated_vapor_pressure) * 100)


def print_law() -> str:
    return print_expression(law)


@validate_input(water_vapor_pressure_=water_vapor_pressure,
    saturated_vapor_pressure_=saturated_vapor_pressure)
@validate_output(relative_humidity)
def calculate_relative_humidity(water_vapor_pressure_: Quantity,
    saturated_vapor_pressure_: Quantity) -> float:
    result_expr = solve(law, relative_humidity, dict=True)[0][relative_humidity]
    result_expr = result_expr.subs({
        water_vapor_pressure: water_vapor_pressure_,
        saturated_vapor_pressure: saturated_vapor_pressure_
    })
    return float(convert_to(Quantity(result_expr), S.One).evalf())
