from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
# Kinetic energy of body EK = (m * v^2) / 2
# where:
# m - body mass
# v - body velocity
kinetic_energy_of_body, body_mass, body_velocity = symbols('kinetic_energy_of_body body_mass body_velocity')
law = Eq(kinetic_energy_of_body, body_mass * body_velocity**2 / 2)

def print():
    return pretty(law, use_unicode=False)

@validate_input(body_mass_=units.mass, body_velocity_=units.velocity)
@validate_output(units.energy)
def calculate_kinetic_energy(body_mass_: Quantity, body_velocity_: Quantity) -> Quantity:
    result_energy_expr = solve(law, kinetic_energy_of_body, dict=True)[0][kinetic_energy_of_body]
    result_expr = result_energy_expr.subs({body_mass: body_mass_, body_velocity: body_velocity_})
    return expr_to_quantity(result_expr, 'kinetic_energy_of_body')