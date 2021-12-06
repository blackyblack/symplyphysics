from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)
from symplyphysics.laws.thermodynamics import pressure_from_temperature_and_volume as thermodynamics_law

# Description
## Boyle's law (Isothermal process): T = const, P1 * V1 = P2 * V2
## Where:
## P is pressure,
## V is volume,
## T is temperature
## V1 is initial volume, V2 is resulting volume
## P1 is initial pressure, P1 is resulting pressure

pressure_start, pressure_end = symbols('pressure_start pressure_end')
volume_start, volume_end = symbols('volume_start volume_end')
law = Eq(pressure_start * volume_start, pressure_end * volume_end)

## Derive the same law from the general ideal gas law

temperature_start, temperature_end = symbols('temperature_start temperature_end')

isothermal_condition = Eq(temperature_start, temperature_end)

eq_start = thermodynamics_law.law.subs({
    thermodynamics_law.temperature: temperature_start,
    thermodynamics_law.volume: volume_start,
    thermodynamics_law.pressure: pressure_start})

eq_end = thermodynamics_law.law.subs({
    thermodynamics_law.temperature: temperature_end,
    thermodynamics_law.volume: volume_end,
    thermodynamics_law.pressure: pressure_end})

derived_law = [eq_start, eq_end, isothermal_condition]

## Check the equivalence of 'law' and 'derived_law'
derived_pressure_end = solve(derived_law,
    (temperature_start, temperature_end, pressure_end), dict=True)[0][pressure_end]
assert solve(law, pressure_end, dict=True)[0][pressure_end] == derived_pressure_end

def print():
    return pretty(law, use_unicode=False)

@validate_input(pressure_start_=units.pressure, pressure_end_=units.pressure, volume_start_=units.volume)
@validate_output(units.volume)
def calculate_volume(pressure_start_: Quantity, volume_start_: Quantity, pressure_end_: Quantity) -> Quantity:
    solved = solve(law, volume_end, dict=True)[0][volume_end]
    result_expr = solved.subs({
        pressure_start: pressure_start_,
        volume_start: volume_start_,
        pressure_end: pressure_end_})
    return expr_to_quantity(result_expr, 'volume_end')
