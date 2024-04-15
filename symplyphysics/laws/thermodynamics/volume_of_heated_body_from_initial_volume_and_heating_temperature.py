from sympy import Eq, solve, dsolve, Symbol as SymSymbol
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import volumetric_coefficient_of_thermal_expansion as coef_def

# Description
## As the temperature increases, the intensity of the thermal motion of the particles of matter increases. This leads to the fact that the molecules are more "actively" repelled from each other.
## As the empty spaces between molecules increase, most bodies increase in size when heated.

## Law: V = V_start * (1 + alpha_V * (t_end - t_start))
## Where:
## V_start is body volume at the initial temperature
## V is volume of body at temperature t
## alpha_V is the coefficient of volumetric expansion of the body
## t_start is initial temperature
## t_end is final temperature

## Conditions
## The pressure must be constant.
## There is no phase transition, eg melting or evaporation.
## The temperature change (difference between `t_end` and `t_start`) is small.

final_volume = Symbol("final_volume", units.volume)
start_volume = Symbol("start_volume", units.volume)
expansion_coefficient = Symbol("coefficient", 1 / units.temperature)
start_temperature = Symbol("start_temperature", units.temperature)
final_temperature = Symbol("final_temperature", units.temperature)

law = Eq(final_volume,
    start_volume * (1 + expansion_coefficient * (final_temperature - start_temperature)))

# Derive from [definition of expansion coefficient](../../definitions/volumetric_coefficient_of_thermal_expansion.py)

_volume = coef_def.volume
_temperature = coef_def.temperature
_temperature_change = SymSymbol("temperature_change")

_diff_eqn = coef_def.definition.subs(coef_def.volumetric_expansion_coefficient, expansion_coefficient)

_volume_dsolved = dsolve(_diff_eqn, _volume(_temperature), ics={_volume(start_temperature): start_volume}).rhs

_volume_subs = _volume_dsolved.subs(_temperature, start_temperature + _temperature_change).simplify()

_volume_series = _volume_subs.series(_temperature_change, 0, 2).removeO()

_volume_derived = _volume_series.subs(_temperature_change, final_temperature - start_temperature)

assert expr_equals(_volume_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(start_volume_=start_volume,
    expansion_coefficient_=expansion_coefficient,
    final_temperature_=final_temperature,
    start_temperature_=start_temperature)
@validate_output(final_volume)
def calculate_final_volume(start_volume_: Quantity, expansion_coefficient_: Quantity,
    final_temperature_: Quantity, start_temperature_: Quantity) -> Quantity:
    result_expr = solve(law, final_volume, dict=True)[0][final_volume]
    result_volume = result_expr.subs({
        start_volume: start_volume_,
        expansion_coefficient: expansion_coefficient_,
        final_temperature: final_temperature_,
        start_temperature: start_temperature_
    })
    return Quantity(result_volume)
