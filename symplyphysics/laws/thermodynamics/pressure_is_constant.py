from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)
from . import pressure_from_temperature_and_volume as thermodynamics_law

# Description
## Gay-Lussac's law (Isobaric process): P = const, T1 / V1 = T2 / V2
## Where:
## P is pressure,
## V is volume,
## T is temperature
## V1 is initial volume, V2 is resulting volume
## T1 is initial temperature, T2 is resulting temperature

temperature_start, temperature_end = symbols('temperature_start temperature_end')
volume_start, volume_end = symbols('volume_start volume_end')
law = Eq(temperature_start / volume_start, temperature_end / volume_end)

## Derive the same law from the general ideal gas law
eq_start = thermodynamics_law.law.subs({
    thermodynamics_law.temperature: temperature_start,
    thermodynamics_law.volume: volume_start})
eq_end = thermodynamics_law.law.subs({
    thermodynamics_law.temperature: temperature_end,
    thermodynamics_law.volume: volume_end})
temperature_start_law = solve([eq_start, eq_end], (thermodynamics_law.pressure, temperature_start))[temperature_start]
derived_law = Eq(temperature_start, temperature_start_law)

## Check the equivalence of 'law' and 'derived_law'
assert solve(law, temperature_end) == solve(derived_law, temperature_end)

def print():
    return pretty(law, use_unicode=False)

@validate_input(temperature_start_=units.temperature, temperature_end_=units.temperature, volume_start_=units.volume)
@validate_output(units.volume)
def calculate_volume(temperature_start_: Quantity, volume_start_: Quantity, temperature_end_: Quantity) -> Quantity:
    result_expr = solve(law.subs({
        temperature_start: temperature_start_,
        volume_start: volume_start_,
        temperature_end: temperature_end_}))[0]
    return expr_to_quantity(result_expr, 'volume_end')
