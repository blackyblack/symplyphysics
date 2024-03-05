from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless)

# Description
## Relative humidity is the ratio of water vapor pressure p to saturated vapor pressure p0 at a given temperature.
## Saturated vapor is vapor in dynamic equilibrium with a liquid or solid of the same composition in a closed system.
## In such a pair, the number of evaporating molecules is equal to the number of condensing molecules per unit of time.
## Water vapor is the gaseous aggregate state of water. This is the vapor for which humidity is measured.
## Pressures on both vapors should be measured the same way, at the same temperature.
## There is also an oversaturated vapor. Supersaturated vapor is vapor whose pressure exceeds the pressure of saturated
## vapor at a given temperature. It can be obtained by increasing the vapor pressure in a volume free from condensation
## centers (dust particles, ions, small liquid droplets, etc.). Another method of production is cooling saturated vapor
## under the same conditions.

## Law is: ph = p / p_sat, where
## ph - relative humidity,
## p - water vapor static pressure (it is called actual vapor pressure),
## p_sat - saturated vapor pressure.

relative_humidity = Symbol("relative_humidity", dimensionless)

water_vapor_pressure = Symbol("water_vapor_pressure", units.pressure)
saturated_vapor_pressure = Symbol("saturated_vapor_pressure", units.pressure)

law = Eq(relative_humidity, water_vapor_pressure / saturated_vapor_pressure)


def print_law() -> str:
    return print_expression(law)


@validate_input(water_vapor_pressure_=water_vapor_pressure,
    saturated_vapor_pressure_=saturated_vapor_pressure)
@validate_output(relative_humidity)
def calculate_relative_humidity(water_vapor_pressure_: Quantity,
    saturated_vapor_pressure_: Quantity) -> Quantity:
    result_expr = solve(law, relative_humidity, dict=True)[0][relative_humidity]
    result_expr = result_expr.subs({
        water_vapor_pressure: water_vapor_pressure_,
        saturated_vapor_pressure: saturated_vapor_pressure_
    })
    return Quantity(result_expr)
