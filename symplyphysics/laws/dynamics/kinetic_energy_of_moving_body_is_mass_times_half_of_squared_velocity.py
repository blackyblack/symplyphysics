from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## If an object has a mass and moves with some velocity, it bears kinetic energy.
## Law: E = m * V**2 / 2, where
## E is kinetic energy of moving object
## m is mass of this object
## V is linear velocity

kinetic_energy, object_mass, linear_velocity = symbols('kinetic_energy object_mass linear_velocity')
law = Eq(kinetic_energy, object_mass * linear_velocity**2 / 2)

def print():
    return pretty(law, use_unicode=False)

@validate_input(mass_=units.mass, velocity_=units.velocity)
@validate_output(units.energy)
def calculate_energy(mass_: Quantity, velocity_: Quantity) -> Quantity:
    result_energy_expr = solve(law, kinetic_energy, dict=True)[0][kinetic_energy]
    result_expr = result_energy_expr.subs({object_mass: mass_, linear_velocity: velocity_})
    return expr_to_quantity(result_expr, 'energy')
