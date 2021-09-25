from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, assert_equivalent_dimension, expr_to_quantity
)
from . import pressure_from_temperature_and_volume as thermodynamics_law

# Description
## Adiabatic process: P * V^y = const
## Where:
## P is pressure,
## V is volume,
## y is the ratio of specific heats

adiabatic_coefficient = symbols('adiabatic_coefficient')
specific_heats_ratio = symbols('specific_heats_ratio')

adiabatic_condition = Eq(thermodynamics_law.pressure * thermodynamics_law.volume ** specific_heats_ratio, adiabatic_coefficient)
law = [thermodynamics_law.law, adiabatic_condition]

def print():
    return pretty(law, use_unicode=False)

@validate_input(pressure_=units.pressure, volume_=units.volume)
def calculate_adiabatic_coefficient(
  pressure_: Quantity,
  volume_: Quantity,
  specific_heats_ratio_: float) :
    result_coefficient_expr = solve(law,
      (thermodynamics_law.temperature, adiabatic_coefficient))[adiabatic_coefficient]

    result_expr = result_coefficient_expr.subs({
        thermodynamics_law.pressure: pressure_,
        thermodynamics_law.volume: volume_,
        specific_heats_ratio: specific_heats_ratio_})
    result = expr_to_quantity(result_expr, 'adiabatic_coefficient')
    assert_equivalent_dimension(
      result,
      'validate_output',
      'return',
      'calculate_adiabatic_coefficient',
      units.pressure * units.volume**specific_heats_ratio_)
    return result
