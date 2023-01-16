from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity,
)

# Description
## Potential energy of body EP = m * g * h
## Where:
## m - body mass
## h - height from Earth surface
## g - free fall acceleration

free_fall_acceleration = units.acceleration_due_to_gravity
potential_energy_of_body, height, body_mass = symbols('potential_energy_of_body height body_mass')
law = Eq(potential_energy_of_body, body_mass * free_fall_acceleration * height)

def print():
    return pretty(law, use_unicode=False)

@validate_input(body_mass_=units.mass, height_=units.length)
@validate_output(units.energy)
def calculate_potential_energy(body_mass_: Quantity, height_: Quantity) -> Quantity:
    result_energy_expr = solve(law, potential_energy_of_body, dict=True)[0][potential_energy_of_body]
    result_expr = result_energy_expr.subs({body_mass: body_mass_, height: height_})
    return expr_to_quantity(result_expr, 'potential_energy_of_body')