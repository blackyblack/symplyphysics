from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, assert_equivalent_dimension, expr_to_quantity
)
from . import pressure_from_temperature_and_volume as thermodynamics_law

# Description
## Adiabatic process: P * V^y = const
## Where:
## P is pressure,
## V is volume,
## y is the ratio of specific heats

specific_heats_ratio = symbols('specific_heats_ratio')
temperature_start, temperature_end = symbols('temperature_start temperature_end')
volume_start, volume_end = symbols('volume_start volume_end')
pressure_start, pressure_end = symbols('pressure_start pressure_end')

adiabatic_condition = Eq(
  pressure_start * (volume_start ** specific_heats_ratio),
  pressure_end * (volume_end ** specific_heats_ratio))

eq_start = thermodynamics_law.law.subs({
  thermodynamics_law.temperature: temperature_start,
  thermodynamics_law.volume: volume_start,
  thermodynamics_law.pressure: pressure_start})

eq_end = thermodynamics_law.law.subs({
  thermodynamics_law.temperature: temperature_end,
  thermodynamics_law.volume: volume_end,
  thermodynamics_law.pressure: pressure_end})

law = [eq_start, eq_end, adiabatic_condition]

def print():
    return pretty(law, use_unicode=False)

@validate_input(mole_count_=units.amount_of_substance, temperature_start_=units.temperature, volume_start_=units.volume, volume_end_=units.volume)
@validate_output(units.pressure)
def calculate_pressure(
  mole_count_: Quantity,
  temperature_start_: Quantity,
  volume_start_: Quantity,
  volume_end_: Quantity,
  specific_heats_ratio_: float) :

    result_pressure_expr = solve(law, (pressure_start, temperature_end, pressure_end))[pressure_end]
    result_pressure = result_pressure_expr.subs({
      thermodynamics_law.mole_count: mole_count_,
      temperature_start: temperature_start_,
      volume_start: volume_start_,
      volume_end: volume_end_,
      specific_heats_ratio: specific_heats_ratio_
    })
    return expr_to_quantity(result_pressure, 'pressure_end')
