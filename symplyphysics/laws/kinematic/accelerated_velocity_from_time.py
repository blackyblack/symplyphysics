from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Accelerated velocity is time dependent and increases with time if acceleration is co-directed with velocity and decreases if they are counter-directed.
## Initial velocity is velocity at the start of observation, when time is 0.

velocity = symbols('velocity')
time = symbols('time')
acceleration = symbols('acceleration')
initial_velocity = symbols('initial_velocity')

law = Eq(velocity, initial_velocity + acceleration * time)

def print():
    return pretty(law, use_unicode=False)
'''
@validate_input(generator_mass_=units.mass, object_mass_=units.mass, distance_between_mass_centers_=units.length)
@validate_output(units.force)
def calculate_force(generator_mass_: Quantity, object_mass_: Quantity, distance_between_mass_centers_: Quantity) -> Quantity:
    result_force_expr = solve(law, gravity_force, dict=True)[0][gravity_force]
    result_expr = result_force_expr.subs({generator_mass: generator_mass_, object_mass: object_mass_, distance_between_mass_centers: distance_between_mass_centers_})
    return expr_to_quantity(result_expr, 'gravity_force')
'''