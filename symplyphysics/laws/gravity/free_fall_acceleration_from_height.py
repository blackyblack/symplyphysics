from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
# Acceleration free fall law: g(R) = (G * M) / R**2
# where :
# G - universal gravity constant  6.672e-11 N*m^2/kg^2
# M - gravity generator mass
# R - distance from the center of the generator
acceleration_free_fall, generator_mass, generator_radius = symbols('acceleration_free_fall generator_mass generator_radius')
law = Eq(acceleration_free_fall, units.gravitational_constant * generator_mass / generator_radius**2)

def print():
    return pretty(law, use_unicode=False)

@validate_input(generator_mass_=units.mass, generator_radius_=units.length)
@validate_output(units.acceleration)
def calculate_acceleration(generator_mass_: Quantity, generator_radius_: Quantity) -> Quantity:
    result_accel_expr = solve(law, acceleration_free_fall, dict=True)[0][acceleration_free_fall]
    result_expr = result_accel_expr.subs({generator_mass: generator_mass_, generator_radius: generator_radius_})
    return expr_to_quantity(result_expr, 'acceleration_free_fall')