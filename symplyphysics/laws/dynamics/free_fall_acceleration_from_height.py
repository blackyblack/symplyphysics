from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
# Accelerate free fall law: g(h)=(G*M)/(R+h)**2
# where :
# G - universal gravity constant  6.672e-11 N*m^2/kg^2
# M - Earth mass constant         5.976e+24 kg
# R - Earth radius constant       6.371e+6 m
# h - height

acceleration_free_fall, constant_gravitation, earth_mass, earth_radius, height = symbols('acceleration_free_fall constant_gravitation earth_mass earth_radius height')
law = Eq(acceleration_free_fall, constant_gravitation*earth_mass / (earth_radius+height)**2)

def print():
    return pretty(law, use_unicode=False)

@validate_input(height_=units.meter)
@validate_output(units.acceleration)
def calculate_acceleration(height_: Quantity) -> Quantity:
    result_accel_expr = solve(law, acceleration_free_fall, dict=True)[0][acceleration_free_fall]
    result_expr = result_accel_expr.subs({height: height_})
    return expr_to_quantity(result_expr, 'acceleration_free_fall')
