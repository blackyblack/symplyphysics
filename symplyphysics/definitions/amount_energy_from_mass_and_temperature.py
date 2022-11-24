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
amount_energy, specific_heat_capacity, body_mass, final_temperature, initial_temperature = symbols(
    'amount_energy specific_heat_capacity body_mass final_temperature initial_temperature'
)
law = Eq(amount_energy, specific_heat_capacity * body_mass * (final_temperature - initial_temperature))

def print():
    return pretty(law, use_unicode=False)

@validate_input(specific_heat_capacity_=units.energy / (units.mass * units.temperature), body_mass_=units.mass,
    final_temperature_=units.temperature, initial_temperature_=units.temperature)
@validate_output(units.energy)
def calculate_amount_energy(specific_heat_capacity_: Quantity, body_mass_: Quantity, final_temperature_: Quantity,
    initial_temperature_: Quantity) -> Quantity:
    result_amount_energy_expr = solve(law, amount_energy_for_body_heat, dict=True)[0][amount_energy_for_body_heat]
    result_expr = result_amount_energy_expr.subs({specific_heat_capacity: specific_heat_capacity_,
                body_mass: body_mass_, final_temperature: final_temperature_, initial_temperature: initial_temperature_})
    return expr_to_quantity(result_expr, 'amount_energy')