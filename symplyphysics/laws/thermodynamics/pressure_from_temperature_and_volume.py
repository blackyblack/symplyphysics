from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Ideal gas law: P * V = n * R * T
## Where:
## P is pressure,
## V is volume,
## n is number of moles,
## R is ideal gas constant,
## T is temperature

pressure, volume, mole_count, temperature = symbols('pressure volume mole_count temperature')
law = Eq(pressure, mole_count * temperature * units.molar_gas_constant / volume)

def print():
    return pretty(law, use_unicode=False)

@validate_input(volume_=units.volume, temperature_=units.temperature, mole_count_=units.amount_of_substance)
@validate_output(units.pressure)
def calculate_pressure(volume_: Quantity, temperature_: Quantity, mole_count_: Quantity) -> Quantity:
    result_pressure_expr = solve(law, pressure)[0]
    result_expr = result_pressure_expr.subs({
        volume: volume_,
        temperature: temperature_,
        mole_count: mole_count_})
    return expr_to_quantity(result_expr, 'pressure')