from sympy import Eq, solve
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    quantities,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import volume_is_constant as isochoric_law

# Description
## The change in gas pressure when the temperature changes relative to 273 kelvin depends on temperature,
## standard pressure and thermal coefficient.
## https://www.fxyz.ru/формулы_по_физике/термодинамика_теория_теплоты/тепловое_расширение/расширение_газа/термический_коэффициент_давления/

## Law is: dp = p0 * (y * T - 1), where
## dp - pressure change,
## p0 - pressure at zero degrees celsius,
## y - thermal coefficient,
## T - temperature.

# Conditions:
## - the volume of gas remains constant when heated

pressure_change = Symbol("pressure", units.pressure)

standard_pressure = Symbol("standard_pressure", units.pressure)
thermal_coefficient = 1 / quantities.standard_conditions_temperature
temperature = symbols.thermodynamics.temperature

law = Eq(pressure_change, standard_pressure * (thermal_coefficient * temperature - 1))

# Derive law from Charles' law of isochoric heating of gas.

_isochoric_eqn = isochoric_law.law.subs({
    isochoric_law.pressure_start: standard_pressure,
    isochoric_law.pressure_end: standard_pressure + pressure_change,
    isochoric_law.temperature_start: quantities.standard_conditions_temperature,
    isochoric_law.temperature_end: temperature,
})

_pressure_change_expr = solve(_isochoric_eqn, pressure_change)[0]

assert expr_equals(_pressure_change_expr, law.rhs)


@validate_input(standard_pressure_=standard_pressure, temperature_=temperature)
@validate_output(pressure_change)
def calculate_pressure_change(standard_pressure_: Quantity, temperature_: Quantity) -> Quantity:
    result_expr = law.rhs.subs({standard_pressure: standard_pressure_, temperature: temperature_})
    return Quantity(result_expr)
