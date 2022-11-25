from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
# Amount of energy for body heat Q = C * m * (t2 - t1)
# where:
# C - specific heat capacity
# m - body mass
# t2 - final temperature
# t1 - initial temperature
# Note: specific heat capacity C is parameter indicates how much energy in joules
# must be expended to heat a unit mass of a substance by one unit temperature.
amount_energy, specific_heat_capacity, body_mass, temperature_begin, temperature_end = symbols(
    'amount_energy specific_heat_capacity body_mass temperature_begin temperature_end'
)
law = Eq(amount_energy, specific_heat_capacity * body_mass * (temperature_begin - temperature_end))

def print():
    return pretty(law, use_unicode=False)

@validate_input(specific_heat_capacity_=units.energy / (units.mass * units.temperature), body_mass_=units.mass,
    temperature_begin_=units.temperature, temperature_end_=units.temperature,
)
@validate_output(units.energy)
def calculate_amount_energy(specific_heat_capacity_: Quantity, body_mass_: Quantity, temperature_begin_: Quantity,
    temperature_end_: Quantity) -> Quantity:
    result_amount_energy_expr = solve(law, amount_energy, dict=True)[0][amount_energy]
    result_expr = result_amount_energy_expr.subs({specific_heat_capacity: specific_heat_capacity_,
                body_mass: body_mass_, temperature_begin: temperature_begin_, temperature_end: temperature_end_})
    return expr_to_quantity(result_expr, 'amount_energy')