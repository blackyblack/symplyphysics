from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)
from symplyphysics.laws.thermodynamics import pressure_from_temperature_and_volume as thermodynamics_law

# Description
## Isochoric process: V = const, P1 * T2 = P2 * T1
## Where:
## P is pressure,
## V is volume,
## T is temperature
## T1 is initial temperature, T2 is resulting temperature
## P1 is initial pressure, P1 is resulting pressure

pressure_start, pressure_end = symbols('pressure_start pressure_end')
temperature_start, temperature_end = symbols('temperature_start temperature_end')
law = Eq(pressure_start * temperature_end, pressure_end * temperature_start)

## Derive the same law from the general ideal gas law

volume_start, volume_end = symbols('volume_start volume_end')

isochoric_condition = Eq(volume_start, volume_end)

eq_start = thermodynamics_law.law.subs({
    thermodynamics_law.temperature: temperature_start,
    thermodynamics_law.volume: volume_start,
    thermodynamics_law.pressure: pressure_start})

eq_end = thermodynamics_law.law.subs({
    thermodynamics_law.temperature: temperature_end,
    thermodynamics_law.volume: volume_end,
    thermodynamics_law.pressure: pressure_end})

derived_law = [eq_start, eq_end, isochoric_condition]

## Check the equivalence of 'law' and 'derived_law'
derived_pressure_end = solve(derived_law, (volume_start, volume_end, pressure_end), dict=True)[0][pressure_end]
assert solve(law, pressure_end, dict=True)[0][pressure_end] == derived_pressure_end

def print():
    return pretty(law, use_unicode=False)

@validate_input(temperature_start_=units.temperature, pressure_start_=units.pressure, temperature_end_=units.temperature)
@validate_output(units.pressure)
def calculate_pressure(temperature_start_: Quantity, pressure_start_: Quantity, temperature_end_: Quantity) -> Quantity:
    solved = solve(law, pressure_end, dict=True)[0][pressure_end]
    result_expr = solved.subs({
        pressure_start: pressure_start_,
        temperature_start: temperature_start_,
        temperature_end: temperature_end_})
    return expr_to_quantity(result_expr, 'pressure_end')
