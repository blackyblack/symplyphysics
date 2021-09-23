from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)
from . import pressure_from_temperature_and_volume as thermodynamics_law

# Description
## Boyle's law (Isotermic process): T = const, P1 * V1 = P2 * V2
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
eq_start = thermodynamics_law.law.subs({
    thermodynamics_law.pressure: pressure_start,
    thermodynamics_law.volume: volume_start})
eq_end = thermodynamics_law.law.subs({
    thermodynamics_law.pressure: pressure_end,
    thermodynamics_law.volume: volume_end})
pressure_start_law = solve([eq_start, eq_end], (thermodynamics_law.temperature, pressure_start))[pressure_start]
derived_law = Eq(pressure_start, pressure_start_law)

## Check the equivalence of 'law' and 'derived_law'
assert solve(law, pressure_end) == solve(derived_law, pressure_end)

def print():
    return pretty(law, use_unicode=False)

@validate_input(pressure_start_=units.pressure, pressure_end_=units.pressure, volume_start_=units.volume)
@validate_output(units.volume)
def calculate_volume(pressure_start_: Quantity, volume_start_: Quantity, pressure_end_: Quantity) -> Quantity:
    result_expr = solve(law.subs({
        pressure_start: pressure_start_,
        volume_start: volume_start_,
        pressure_end: pressure_end_}))[0]
    return expr_to_quantity(result_expr, 'volume_end')
