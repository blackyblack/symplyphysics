from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity, SI, sqrt
)

# Description
## If hole appears in side wall or bottom of tank with liquid, liquid starts flowing out of this tank with some velocity. 
## Velocity might be alculated with help of Torrichelli's formula. 
## Law: v = sqrt(2 * g * h), where
## v is velocity of liquid flowing out of the hole
## g is gravitation acceleration
## h is height of liquid above the hole.

# Conditions
## Liquid is ideal. Hole is small.

liquid_velocity, gravitation_acceleration, height_above_hole = symbols('liquid_velocity gravitation_acceleration height_above_hole')
law = Eq(liquid_velocity, sqrt(2 * gravitation_acceleration * height_above_hole))

def print():
    return pretty(law, use_unicode=False)

@validate_input(gravitation_acc_ = units.length / units.time**2, height_ = units.length)
@validate_output(units.velocity)
def calculate_velocity(gravitation_acc_: Quantity, height_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, liquid_velocity, dict=True)[0][liquid_velocity]        
    result_expr = result_velocity_expr.subs({gravitation_acceleration: gravitation_acc_, height_above_hole: height_})
    return expr_to_quantity(result_expr, 'liquid_energy')
